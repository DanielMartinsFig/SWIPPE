library(shiny)
library(bslib)
library(dplyr)
library(ggplot2)
library(scales)
library(purrr)
library(viridis)

# ------------------ Helpers for SWIPPE (same as before) ------------------
sigma_y_turner <- function(X, stab) {
  ifelse(stab == 1, 0.22*X*(1+0.0001*X)^-0.5,
         ifelse(stab == 2, 0.16*X*(1+0.0001*X)^-0.5,
                ifelse(stab == 3, 0.11*X*(1+0.0001*X)^-0.5,
                       ifelse(stab == 4, 0.08*X*(1+0.0001*X)^-0.5,
                              ifelse(stab == 5, 0.06*X*(1+0.0001*X)^-0.5,
                                     0.04*X*(1+0.0001*X)^-0.5)))))
}
sigma_z_turner <- function(X, stab) {
  ifelse(stab == 1, 0.20*X,
         ifelse(stab == 2, 0.12*X,
                ifelse(stab == 3, 0.08*X*(1+0.0002*X)^-0.5,
                       ifelse(stab == 4, 0.06*X*(1+0.0015*X)^-0.5,
                              ifelse(stab == 5, 0.03*X*(1+0.0003*X)^-1,
                                     0.016*X*(1+0.0003*X)^-1)))))
}
clamp01 <- function(x) pmax(0, pmin(1, x))
norm_texture <- function(sand, silt, clay) {
  s <- sand + silt + clay
  if (is.na(s) || s == 0) return(c(sand = 0, silt = 0, clay = 0))
  c(sand = sand/s, silt = silt/s, clay = clay/s)
}

compute_once <- function(
    DU10, DIAM, sand, silt, clay, Foc, ST, TS,
    RRUF, wn2, RHTT, Fconventional, Fconservation, Fnotill,
    Fresidue, Fcover, Ccrop, THW, ANG, FL, FW,
    Koc, m, Stab, X, Z, Hs, reflect_on = TRUE
){
  tex <- norm_texture(sand, silt, clay) * 100
  sand <- tex["sand"]; silt <- tex["silt"]; clay <- tex["clay"]
  
  Fom <- Foc * 1.724
  USTR  <- 0.0408 * DU10
  USTRT <- 0.0161 * sqrt(DIAM)
  
  WPt <- -0.024*sand + 0.487*clay + 0.006*Fom + 0.005*(sand*Fom) -
    -0.013*(clay*Fom) + 0.068*(sand*clay) + 0.031
  WP <- pmax(WPt + ((0.14*WPt) - 0.02), 1e-6)
  
  A <- pmax(0, USTR^2 - USTRT^2 - 0.5*(ST/WP)^2)
  YWR <- TS * 0.255 * (A)^1.5
  
  SEF <- 29.09 + 0.31*sand + 0.17*silt + 0.33*(sand/clay) - 2.59*Fom
  RRF <- 11.9 * (1 - exp(-((RRUF/9.8)^1.3)))
  RIF <- abs(sin(wn2*pi/180)) * 1.27 * (RHTT^0.52)
  RFB <- RRF + RIF
  RFC <- 0.77 * (1.002^RHTT)
  FRF <- 1 - exp(-(((10*(pi/180))/pmax(RFB, 1e-6))^RFC))
  
  Fconventional <- clamp01(Fconventional); Fconservation <- clamp01(Fconservation); Fnotill <- clamp01(Fnotill)
  Fresidue <- clamp01(Fresidue); Fcover <- clamp01(Fcover)
  Ctillage <- Fconventional + Fconservation*0.35 + Fnotill*0.25
  Cresidue <- (Fresidue*0.88) + (1 - Fresidue)
  Ccover   <- Fcover*0.8 + (1 - Fcover)
  Cmanagement <- Ctillage*Cresidue*Ccover
  VGF <- Ccrop * Cmanagement
  
  BT <- 1.571 + (THW*pi/180) - (ANG*pi/180)
  WL <- FL*FW / (FL*abs(cos(BT)) + FW*abs(sin(BT)))
  WL <- pmax(WL, 1e-6)
  FD <- 1 - exp(-((WL/0.07)^2))
  
  YW <- SEF * FRF * VGF * FD * (YWR / WL)
  
  Kd <- Koc * Foc
  Vw <- 460 * ST
  denom <- pmax(460*Kd + Vw, 1e-12)
  Cps <- (m * Kd) / denom            # mg/kg
  Qdust <- YW / TS                   # kg/s per m2
  Qpesticide <- Qdust * Cps * 1e-6   # kg/s (mg->kg)
  
  # Plume centerline at (X,Z)
  U <- ifelse(DU10 < 1, NA_real_, DU10)
  Vp <- sqrt((4/3)*((0.00016*9.8*(1000-1.225))/(12*1.225)))
  sy <- sigma_y_turner(X, Stab)
  sz <- sigma_z_turner(X, Stab)
  Heff <- Hs - (Vp*X)/U
  vert_core <- exp(-((Z - Heff)^2)/(2*sz^2))
  vert_term <- if (isTRUE(reflect_on)) vert_core + exp(-((Z + Heff)^2)/(2*sz^2)) else vert_core
  base <- (vert_term) / (2*pi*sy*sz*U)
  Cdust <- Qdust * base
  Cpesticide <- Qpesticide * base
  
  tibble(
    YW = YW, Qdust = Qdust, Cps_mg_per_kg = Cps,
    Cdust_ng_m3 = Cdust * 1e12, Cpesticide_ng_m3 = Cpesticide * 1e12,
    WL_km = WL, FRF = FRF, VGF = VGF, FD = FD, SEF = SEF,
    U_valid = !is.na(U)
  )
}

# ---------------- Temporal plume (your PG moving-frame version) --------------
# Pasquill–Gifford rural parameters
pg_params <- list(
  A = list(a=0.22, b=0.894, c=0.20,  d=1.0),
  B = list(a=0.16, b=0.894, c=0.12,  d=1.0),
  C = list(a=0.11, b=0.894, c=0.08,  d=0.915),
  D = list(a=0.08, b=0.894, c=0.06,  d=0.914),
  E = list(a=0.06, b=0.894, c=0.03,  d=0.899),
  F = list(a=0.04, b=0.894, c=0.016, d=0.902)
)
pg_sigmas <- function(x, cls) {
  p <- pg_params[[cls]]
  x <- pmax(x, 1e-6)
  sy <- p$a * x^p$b
  sz <- p$c * x^p$d
  list(sy=sy, sz=sz)
}

# Gaussian plume (ground reflection) in moving frame; H_eff option
plume_conc_pg <- function(Q, U, H_eff, x_rel, y, cls) {
  C <- rep(NA_real_, length(x_rel))
  pos <- which(x_rel > 0)
  if (!length(pos)) return(C)
  sg <- pg_sigmas(x_rel[pos], cls)
  sy <- sg$sy; sz <- sg$sz
  core <- Q / (2*pi*U*sy*sz)
  cy <- exp(-(y[pos]^2)/(2*sy^2))
  cz <- 2 * exp(-(H_eff^2)/(2*sz^2))  # image source
  C[pos] <- core * cy * cz
  C
}

# ----------------------------- UI -------------------------------------------
theme <- bs_theme(bootswatch = "flatly")
ui <- page_fluid(
  theme = theme,
  tags$head(
    tags$style(HTML("
      .control-card { min-width: 280px; border: none; }
      .control-card .card-body { padding: 12px 14px; }
      .control-card h4 { margin-top: 4px; font-size: 1.05rem; }
      .container-fluid { max-width: 1400px; }
      .bg-met   { background-color: #e6f2ff; }
      .bg-soil  { background-color: #f4e7d3; }
      .bg-veg   { background-color: #e6f7e8; }
      .bg-geom  { background-color: #efe6ff; }
      .control-card .card-body > * { margin-bottom: 6px; }
      .small-note { font-size: 12px; opacity: 0.8; }
      .outputs { margin-top: 18px; }
    "))
  ),
  
  titlePanel("SWIPPE — Wind erosion of particle-phase pesticides and temporal plume visualizer)"),
  
  navset_tab(
    nav_panel("Calculator",
              div(
                # --- New instruction box ---
                card(
                  class = "mb-3 border-info",
                  card_body(
                    tags$div(
                      style = "background-color:#e8f4fd; border-left:4px solid #3399ff; padding:10px 14px; border-radius:6px;",
                      HTML("<b>Instructions:</b> Fill in and choose the different values for each variable, 
                        then click the <b>Run scenario</b> button below to compute results.")
                    )
                  )
                ),
                
                # --- Original heading ---
                h3("Controls"),
                layout_columns(col_widths = c(3,3,3,3),
                               card(class = "control-card bg-met",
                                    card_header("Meteorology"),
                                    card_body(
                                      sliderInput("DU10", "Wind @10 m (m/s)", min = 0, max = 20, value = 5, step = 0.1),
                                      selectInput("Stab", "Stability class (Turner/steady)", setNames(1:6, c("A","B","C","D","E","F")), selected = 4),
                                      numericInput("THW", "Wind dir from N (°)", 0, 0, 360, 1),
                                      numericInput("ANG", "Field length angle from N (°)", 0, 0, 360, 1),
                                      checkboxInput("reflect_on", "Include ground reflection", TRUE)
                                    )
                               ),
                               card(class = "control-card bg-soil",
                                    card_header("Soil & texture (%)"),
                                    card_body(
                                      numericInput("sand", "Sand (%)", 60, 0, 100, 1),
                                      numericInput("silt", "Silt (%)", 25, 0, 100, 1),
                                      numericInput("clay", "Clay (%)", 15, 0, 100, 1),
                                      div(class = "small-note", "Texture auto-normalized to 100%. Sum: ", textOutput("textureSum", inline=TRUE)),
                                      sliderInput("Foc", "Organic carbon, FOC (%)", 0, 60, 1.5, 0.1),
                                      sliderInput("ST", "Topsoil water, ST (%)", 0, 20, 1.5, 0.1),
                                      numericInput("TS", "Time step, TS (s)", 3600, 60, step = 60)
                                    )
                               ),
                               card(class = "control-card bg-veg",
                                    card_header("Roughness & Vegetation"),
                                    card_body(
                                      numericInput("RRUF", "Random roughness RRUF (mm)", 10, 0, step = 1),
                                      numericInput("wn2", "Wind vs ridge angle (°)", 90, 0, 180, 1),
                                      numericInput("RHTT", "Ridge height RHTT (mm)", 30, 0, step = 1),
                                      tags$hr(),
                                      sliderInput("Fconventional", "Conventional tillage frac.", 0, 1, 0.5, 0.01),
                                      sliderInput("Fconservation", "Conservation tillage frac.", 0, 1, 0.3, 0.01),
                                      sliderInput("Fnotill", "No-till frac.", 0, 1, 0.2, 0.01),
                                      sliderInput("Fresidue", "Residue-treated frac.", 0, 1, 0.4, 0.01),
                                      sliderInput("Fcover", "Cover crops frac.", 0, 1, 0.3, 0.01),
                                      numericInput("Ccrop", "Crop C-factor (Panagos 2015)", 0.34, 0, step = 0.01)
                                    )
                               ),
                               card(class = "control-card bg-geom",
                                    card_header("Geometry, Sorption & Receptor"),
                                    card_body(
                                      numericInput("FL", "Field length FL (km)", 0.5, 0.001, step = 0.01),
                                      numericInput("FW", "Field width FW (km)", 0.5, 0.001, step = 0.01),
                                      tags$hr(),
                                      numericInput("Koc", "Koc (mL/g)", 500, 0, step = 10),
                                      numericInput("m", "Pesticide mass per m² (mg)", 5, 0, step = 0.1),
                                      tags$hr(),
                                      numericInput("DIAM", "Soil particle diameter (m)", 75e-6, 1e-6, step = 1e-6),
                                      numericInput("Hs", "Source height Hs (m)", 0, 0, step = 0.1),
                                      numericInput("X", "Receptor distance X (m)", 50, 1, step = 1),
                                      numericInput("Z", "Sampler height Z (m)", 1.5, 0, step = 0.1)
                                    ),
                                    card_footer(actionButton("run1", "Run scenario", class = "btn-primary w-100"))
                               )
                )
              ),
              
              # ---- Results (steady) ----
              div(class = "outputs",
                  h3("Results"),
                  fluidRow(
                    column(4, card(card_header("Dust concentration (ng m⁻³)"), card_body(h2(textOutput("Cdust_txt"))))),
                    column(4, card(card_header("Pesticide in soil (Cps, mg kg⁻¹)"), card_body(h2(textOutput("Cps_txt"))))),
                    column(4, card(card_header("Particle-phase pesticide (ng m⁻³)"), card_body(h2(textOutput("Cpest_txt")))))
                  ),
                  fluidRow(
                    column(6, card(card_header("Plume centerline vs distance (steady)"),
                                   card_body(plotOutput("plumePlot", height = 320)))),
                    column(6, card(card_header("Output"),
                                   card_body(tableOutput("diagTable"))))
                  ),
                  downloadButton("dl_single", "Download single-scenario CSV"),
                  
                  # ---- Temporal plume (PG moving frame) ----
                  tags$hr(),
                  h3("Temporal plume (moving frame: x_rel = x − U·t)"),
                  fluidRow(
                    column(4,
                           selectInput("pg_cls", "PG stability class (temporal)", choices = c("A","B","C","D","E","F"), selected = "D"),
                           sliderInput("pg_U", "Wind speed U (m/s)", 0.5, 15, value = 3, step = 0.1),
                           sliderInput("xmax_km", "Downwind extent (km)", 0.2, 10, value = 3, step = 0.1),
                           sliderInput("ymax_km", "Crosswind half-width (km)", 0.1, 5, value = 1, step = 0.1),
                           sliderInput("res_m", "Grid resolution (m)", 10, 200, value = 50, step = 10),
                           checkboxInput("pg_log", "Log color scale", TRUE),
                           tags$hr(),
                           sliderInput("dt_s", "Time step Δt (s)", 1, 120, value = 10, step = 1),
                           sliderInput("t_window_min", "Simulated time window (min)", min = 0, max = 120, value = c(0, 20), step = 1),
                           sliderInput("t_min", "Time t (min) — use Play ▶", min = 0, max = 60, value = 0, step = 0.5),
                           sliderInput("focusX", "Focus distance X (m)", 5, 500, value = 50, step = 5),
                           sliderInput("play_ms", "Play speed (ms)", 100, 1500, value = 300, step = 50),
                           fluidRow(
                             column(6, actionButton("btnStart", "Start")),
                             column(6, actionButton("btnStop",  "Stop"))
                           ),
                           helpText("Temporal plume uses your PG spec and advects the steady plume with x_rel = x − U·t.")
                    ),
                    column(8,
                           plotOutput("pg_plan", height = "520px"),
                           plotOutput("pg_centerline", height = "160px"),
                           plotOutput("pg_time_focus", height = "160px")
                    )
                  )
              )
    ),
    
    nav_panel("Inputs explanation",
              div(class = "container-fluid",
                  h3("What each input means (and typical ranges)"),
                  
                  h4("Meteorology"),
                  tags$ul(
                    tags$li(HTML("<b>Wind @10 m (m/s)</b> — Mean wind speed measured at 10 m. Typical near-surface winds range ~0–15 m/s; strong events >15 m/s are possible in storms. Used to compute friction velocity and erosion initiation.")),
                    tags$li(HTML("<b>Stability class (A–F)</b> — Pasquill–Gifford/Turner stability: A,B = very/ moderately unstable (sunny, convective); C = slightly unstable; D = neutral (overcast, windy); E,F = stable (night, light winds). Controls σ<sub>y</sub>, σ<sub>z</sub> growth and reflection term. <span class='small-note'>(Turner, 1969)</span>")),
                    tags$li(HTML("<b>Wind dir from N (°)</b> — Meteorological direction the wind blows from (0–360). Used with field orientation to compute effective fetch.")),
                    tags$li(HTML("<b>Include ground reflection</b> — Image-source term that doubles the vertical Gaussian about z=0 for a ground-level receptor (standard in simple Gaussian plumes). <span class='small-note'>(Turner, 1969)</span>"))
                  ),
                  
                  h4("Soil & texture"),
                  tags$ul(
                    tags$li(HTML("<b>Sand / Silt / Clay (%)</b> — Particle size fractions; auto-normalized to sum 100%. USDA texture triangle is the reference classification. <span class='small-note'>(USDA NRCS, 2017)</span>")),
                    tags$li(HTML("<b>Organic carbon, FOC (%)</b> — Soil organic carbon fraction. Commonly ~0.5–5% in mineral topsoils but broader 0.1–10% occurs depending on land use and climate; drives sorption via K<sub>oc</sub>. <span class='small-note'>(USDA NRCS, 2017)</span>")),
                    tags$li(HTML("<b>Topsoil water, ST (%)</b> — Gravimetric water content. Typical agricultural topsoils range roughly 5–30% by mass; dries down under warm/windy conditions.")),
                    tags$li(HTML("<b>Time step, TS (s)</b> — Integration step for wind erosion flux <i>Y<sub>W</sub></i>. Shorter steps resolve transients; 600–3600 s is a practical range."))
                  ),
                  
                  h4("Roughness & Vegetation / Management"),
                  tags$ul(
                    tags$li(HTML("<b>Random roughness RRUF (mm)</b> — Standard deviation of microrelief heights after removing slope/tillage trend; typical freshly tilled fields ~10–30 mm, smoother crusted surfaces ~2–10 mm. Larger values increase sheltering and reduce erosion. <span class='small-note'>(Potter et&nbsp;al., 1990; Vinci et&nbsp;al., 2020)</span>")),
                    tags$li(HTML("<b>Wind vs ridge angle (°)</b> — Angle between wind and oriented roughness (ridges). Max shelter near 90° (cross-ridge).")),
                    tags$li(HTML("<b>Ridge height RHTT (mm)</b> — Oriented roughness height; higher ridges generally reduce erosion flux via sheltering.")),
                    tags$li(HTML("<b>Conventional / Conservation / No-till fractions</b> — Fractions (0–1) of the field under each tillage system; the model converts to an effective management multiplier.")),
                    tags$li(HTML("<b>Residue-treated fraction</b> — Fraction with residue retained (0–1); residue strongly reduces erodibility and C-factor.")),
                    tags$li(HTML("<b>Cover crops fraction</b> — Fraction with live cover in the erosive season (0–1); reduces C-factor.")),
                    tags$li(HTML("<b>Crop C-factor (Panagos 2015)</b> — RUSLE cover-management factor (unitless). Typical values: perennial/forest ≪ 0.1; small grains/forage ~0.05–0.3; row crops higher depending on residue and cover. You can use Panagos et&nbsp;al. (2015) EU-scale estimates as guidance. <span class='small-note'>(Panagos et&nbsp;al., 2015; ESDAC)</span>"))
                  ),
                  
                  h4("Geometry, Sorption & Receptor"),
                  tags$ul(
                    tags$li(HTML("<b>Field length FL (km)</b>, <b>Field width FW (km)</b> — Planform dimensions; combined with wind/row angle to compute effective fetch width <i>WL</i>. Typical fields: 0.1–1 km scale.")),
                    tags$li(HTML("<b>K<sub>oc</sub> (mL/g)</b> — Organic carbon–water partition coefficient driving pesticide sorption (K<sub>d</sub>=K<sub>oc</sub>·f<sub>oc</sub>). Non-ionic pesticides span orders of magnitude (~10–10<sup>5</sup> mL/g). <span class='small-note'>(ECETOC, 2013; Baker & Hites, 2000)</span>")),
                    tags$li(HTML("<b>Pesticide mass per m², m (mg)</b> — Mass of active ingredient present in the wind-erodible surface layer per ground area (mg/m²). Derive from application rate × interception × degradation/volatilization losses.")),
                    tags$li(HTML("<b>Soil particle diameter DIAM (m)</b> — Characteristic eroding particle/aggregate diameter. Wind-erodible fractions often tens to a few hundred micrometres (e.g., 20–200&nbsp;µm); emitted PM includes finer particles.")),
                    tags$li(HTML("<b>Source height Hs (m)</b> — Effective release height (0–2 m typical for surface erosion; higher if lofted).")),
                    tags$li(HTML("<b>Receptor distance X (m)</b> — Downwind distance to sampler; tens to thousands of meters are common in near-field assessments.")),
                    tags$li(HTML("<b>Sampler height Z (m)</b> — Inlet height; 1–2 m typical for near-ground samplers."))
                  ),
                  
                  h4("Temporal plume (PG moving frame)"),
                  tags$ul(
                    tags$li(HTML("<b>PG stability class</b> — Same A–F notion but using Pasquill–Gifford rural σ fits for time-varying plan view. <span class='small-note'>(Briggs, 1973; Turner, 1969)</span>")),
                    tags$li(HTML("<b>Wind speed U (m/s)</b> — Advection speed for moving-frame plume (x<sub>rel</sub> = x − U·t).")),
                    tags$li(HTML("<b>Downwind extent (km), Crosswind half-width (km)</b> — Domain size for the plan view.")),
                    tags$li(HTML("<b>Grid resolution (m)</b> — Raster cell size; 25–100 m is a good balance.")),
                    tags$li(HTML("<b>Log color scale</b> — Toggle to reveal gradients over several orders of magnitude.")),
                    tags$li(HTML("<b>Δt (s), Time window (min), Play</b> — Controls for stepping and animating the moving plume.")),
                    tags$li(HTML("<b>Focus distance X (m)</b> — Used to extract a time-series at a fixed downwind location."))
                  ),
                  
                  tags$hr(),
                  h4("Notes on typical ranges and interpretation"),
                  tags$ul(
                    tags$li("Ranges above are indicative; local soils, crops, and weather can differ substantially."),
                    tags$li(HTML("C-factor and management fractions interact multiplicatively; small improvements (residue/cover) can yield order-of-magnitude reductions in <i>Y<sub>W</sub></i>.")),
                    tags$li("Stability strongly affects downwind decay; stable nights (E–F) keep plumes narrow and elevated; unstable days (A–C) dilute faster.")
                  ),
                  
                  tags$hr(),
                  h4("References"),
                  tags$ol(
                    tags$li(HTML("Turner, D. B. (1969). <i>Workbook of Atmospheric Dispersion Estimates</i>, U.S. EPA. ")),
                    tags$li(HTML("Briggs, G. A. (1973). Diffusion Estimation for Small Emissions; Pasquill–Gifford σ<sub>y</sub>, σ<sub>z</sub> equations.")),
                    tags$li(HTML("USDA NRCS (2017). <i>Soil Survey Manual</i>, Ch. 3 — USDA texture triangle.")),
                    tags$li(HTML("Potter, K. N., Zobeck, T. M., & Hagan, L. J. (1990s). Microrelief/Random roughness and wind erosion indices.")),
                    tags$li(HTML("Vinci, A. et&nbsp;al. (2020). Comparative evaluation of random roughness indices.")),
                    tags$li(HTML("Panagos, P. et&nbsp;al. (2015). Estimating the soil erosion cover-management (C) factor at the European scale.")),
                    tags$li(HTML("ESDAC (JRC). C-factor datasets and documentation for Europe.")),
                    tags$li(HTML("ECETOC (2013). Organic carbon–water partition (K<sub>oc</sub>) overview; ranges across chemicals. - Technical report 123")),
                    tags$li(HTML("Baker, J., & Hites, R. (2000). Estimating K<sub>oc</sub> for persistent organic pollutants."))
                  ),
                  
                  # tiny footnote so they know where to find these in-app too
                  div(class = "small-note",
                      "Citations are provided for methodology/context; see each paper/dataset for full details.")
              )
    ),
    
    nav_panel("About / assumptions",
              tags$ul(
                tags$li("Texture entries (sand/silt/clay) auto-normalized to 100%."),
                tags$li("Steady block uses Turner σy/σz; temporal block uses Pasquill–Gifford rural exponents."),
                tags$li("Ground reflection included in both (image source)."),
                tags$li("Units: ng/m³ for air concentrations; mg/kg for soil pesticide.")
              )
    )
  )
)

# ----------------------------- Server ----------------------------------------
server <- function(input, output, session) {
  # texture sum
  output$textureSum <- renderText({
    s <- sum(c(input$sand, input$silt, input$clay), na.rm = TRUE)
    sprintf("%.1f%%", s)
  })
  
  # steady scenario (SWIPPE centerline at X, Z)
  res_single <- eventReactive(input$run1, {
    compute_once(
      DU10=input$DU10, DIAM=input$DIAM, sand=input$sand, silt=input$silt, clay=input$clay,
      Foc=input$Foc, ST=input$ST, TS=input$TS, RRUF=input$RRUF, wn2=input$wn2, RHTT=input$RHTT,
      Fconventional=input$Fconventional, Fconservation=input$Fconservation, Fnotill=input$Fnotill,
      Fresidue=input$Fresidue, Fcover=input$Fcover, Ccrop=input$Ccrop, THW=input$THW, ANG=input$ANG,
      FL=input$FL, FW=input$FW, Koc=input$Koc, m=input$m, Stab=as.integer(input$Stab),
      X=input$X, Z=input$Z, Hs=input$Hs, reflect_on=isTRUE(input$reflect_on)
    )
  }, ignoreInit = TRUE)
  
  output$Cdust_txt <- renderText({ r <- res_single(); if (is.null(r)) "—" else if (!r$U_valid) "NA (U<1 m/s)" else format(r$Cdust_ng_m3, digits=3, scientific=TRUE) })
  output$Cps_txt   <- renderText({ r <- res_single(); if (is.null(r)) "—" else format(r$Cps_mg_per_kg, digits=3, scientific=TRUE) })
  output$Cpest_txt <- renderText({ r <- res_single(); if (is.null(r)) "—" else if (!r$U_valid) "NA (U<1 m/s)" else format(r$Cpesticide_ng_m3, digits=3, scientific=TRUE) })
  
  output$diagTable <- renderTable({
    r <- res_single(); if (is.null(r)) return(NULL)
    r |> transmute(
      `U valid`=ifelse(U_valid,"Yes","No"),
      `WL (km)`=round(WL_km,4), `SEF`=round(SEF,2), `FRF`=round(FRF,3),
      `VGF`=round(VGF,3), `FD`=round(FD,3), `YW (kg/m2 in TS)`=signif(YW,3),
      `Qdust (kg/s)`=signif(Qdust,3)
    )
  })
  
  output$plumePlot <- renderPlot({
    r <- res_single()
    if (is.null(r) || !r$U_valid) { plot.new(); title("Plume plot unavailable (U < 1 m/s)"); return() }
    Xseq <- seq(5, 500, by = 5)
    tmp <- map_dfr(Xseq, function(xx) {
      rr <- compute_once(
        DU10=input$DU10, DIAM=input$DIAM, sand=input$sand, silt=input$silt, clay=input$clay,
        Foc=input$Foc, ST=input$ST, TS=input$TS, RRUF=input$RRUF, wn2=input$wn2, RHTT=input$RHTT,
        Fconventional=input$Fconventional, Fconservation=input$Fconservation, Fnotill=input$Fnotill,
        Fresidue=input$Fresidue, Fcover=input$Fcover, Ccrop=input$Ccrop, THW=input$THW, ANG=input$ANG,
        FL=input$FL, FW=input$FW, Koc=input$Koc, m=input$m, Stab=as.integer(input$Stab),
        X=xx, Z=input$Z, Hs=input$Hs, reflect_on=isTRUE(input$reflect_on)
      )
      tibble(X=xx, Cp_ng_m3=rr$Cpesticide_ng_m3)
    })
    ggplot(tmp, aes(X, Cp_ng_m3)) +
      geom_line() +
      scale_y_continuous(labels = label_scientific()) +
      labs(x="Distance downwind X (m)", y="Particle-phase pesticide (ng m⁻³)", title="Plume centerline (steady)") +
      theme_bw()
  })
  
  output$dl_single <- downloadHandler(
    filename = function() "swippe_single.csv",
    content = function(file) { write.csv(res_single(), file, row.names = FALSE) }
  )
  
  # ------------------------ Temporal PG animation ---------------------------
  # keep time slider synced to dt + window
  observe({
    updateSliderInput(session, "t_min",
                      min = input$t_window_min[1], max = input$t_window_min[2],
                      step = input$dt_s/60
    )
  })
  
  playing <- reactiveVal(FALSE)
  observeEvent(input$btnStart, { playing(TRUE) })
  observeEvent(input$btnStop,  { playing(FALSE) })
  
  observe({
    if (!playing()) return()
    isolate({
      tnow <- input$t_min
      tmax <- input$t_window_min[2]
      step_min <- input$dt_s/60
      tnext <- if (tnow >= tmax) input$t_window_min[1] else min(tnow + step_min, tmax)
      updateSliderInput(session, "t_min", value = tnext)
    })
    invalidateLater(input$play_ms)
  })
  
  # Emission rate Q for temporal plume: derive from current SWIPPE inputs (kg/s -> g/s),
  # but use centerline at X=1m only to pull Qpesticide terms; or keep simple control:
  # We'll compute Qpesticide from inputs at X=1 m.
  qp_gps <- reactive({
    rr <- compute_once(
      DU10=input$DU10, DIAM=input$DIAM, sand=input$sand, silt=input$silt, clay=input$clay,
      Foc=input$Foc, ST=input$ST, TS=input$TS, RRUF=input$RRUF, wn2=input$wn2, RHTT=input$RHTT,
      Fconventional=input$Fconventional, Fconservation=input$Fconservation, Fnotill=input$Fnotill,
      Fresidue=input$Fresidue, Fcover=input$Fcover, Ccrop=input$Ccrop, THW=input$THW, ANG=input$ANG,
      FL=input$FL, FW=input$FW, Koc=input$Koc, m=input$m, Stab=as.integer(input$Stab),
      X=1, Z=input$Z, Hs=input$Hs, reflect_on=TRUE
    )
    # back out Qpesticide (kg/s) from underlying pieces: Qdust * Cps * 1e-6
    # We must recompute Qdust and Cps exactly as in compute_once:
    tex <- norm_texture(input$sand, input$silt, input$clay) * 100
    sand <- tex["sand"]; silt <- tex["silt"]; clay <- tex["clay"]
    Fom <- input$Foc * 1.724
    USTR  <- 0.0408 * input$DU10
    USTRT <- 0.0161 * sqrt(input$DIAM)
    WPt <- -0.024*sand + 0.487*clay + 0.006*Fom + 0.005*(sand*Fom) -
      -0.013*(clay*Fom) + 0.068*(sand*clay) + 0.031
    WP <- pmax(WPt + ((0.14*WPt) - 0.02), 1e-6)
    A <- pmax(0, USTR^2 - USTRT^2 - 0.5*(input$ST/WP)^2)
    YWR <- input$TS * 0.255 * (A)^1.5
    SEF <- 29.09 + 0.31*sand + 0.17*silt + 0.33*(sand/clay) - 2.59*Fom
    RRF <- 11.9 * (1 - exp(-((input$RRUF/9.8)^1.3)))
    RIF <- abs(sin(input$wn2*pi/180)) * 1.27 * (input$RHTT^0.52)
    RFB <- RRF + RIF
    RFC <- 0.77 * (1.002^input$RHTT)
    FRF <- 1 - exp(-(((10*(pi/180))/pmax(RFB, 1e-6))^RFC))
    Ctillage <- clamp01(input$Fconventional) + clamp01(input$Fconservation)*0.35 + clamp01(input$Fnotill)*0.25
    Cresidue <- (clamp01(input$Fresidue)*0.88) + (1 - clamp01(input$Fresidue))
    Ccover   <- clamp01(input$Fcover)*0.8 + (1 - clamp01(input$Fcover))
    VGF <- input$Ccrop * (Ctillage*Cresidue*Ccover)
    BT <- 1.571 + (input$THW*pi/180) - (input$ANG*pi/180)
    WL <- pmax(input$FL*input$FW / (input$FL*abs(cos(BT)) + input$FW*abs(sin(BT))), 1e-6)
    FD <- 1 - exp(-((WL/0.07)^2))
    YW <- SEF * FRF * VGF * FD * (YWR / WL)
    Kd <- input$Koc * input$Foc
    Vw <- 460 * input$ST
    Cps <- (input$m * Kd) / pmax(460*Kd + Vw, 1e-12)   # mg/kg
    Qdust <- YW / input$TS                             # kg/s per m2
    Qpesticide <- Qdust * Cps * 1e-6                   # kg/s
    Qpesticide * 1e3                                   # g/s for PG plume
  })
  
  # Effective height for temporal: include settling tilt relative to x_rel
  # Use H_eff at source Hs; PG plume uses cz = 2*exp(-H_eff^2/(2*sz^2)).
  # To mimic tilt, we pass x_rel to sigmas and keep Hs fixed (common simplification).
  # If you want explicit tilt: replace H_eff with Hs - (Vp/U) * x_rel (would require per-grid eval).
  
  # Plan view (raster/contours)
  output$pg_plan <- renderPlot({
    xmax <- input$xmax_km * 1000
    ymax <- input$ymax_km * 1000
    dx   <- input$res_m
    
    t_sec <- input$t_min * 60
    x0 <- input$pg_U * t_sec   # advection distance of plume core
    
    xs <- seq(0, xmax, by = dx)
    ys <- seq(-ymax, ymax, by = dx)
    grid <- expand.grid(x = xs, y = ys)
    
    x_rel <- grid$x - x0
    C_g_m3 <- plume_conc_pg(
      Q = qp_gps(), U = input$pg_U, H_eff = input$Hs,
      x_rel = x_rel, y = grid$y, cls = input$pg_cls
    )
    
    df <- data.frame(x = grid$x/1000, y = grid$y/1000, C = C_g_m3)
    df$C_plot <- df$C
    df$C_plot[df$C_plot <= 0 | is.na(df$C_plot)] <- NA_real_
    
    ggplot(df, aes(x, y, z = if (input$pg_log) log10(C_plot) else C_plot)) +
      geom_raster(aes(fill = if (input$pg_log) log10(C_plot) else C_plot), na.rm = TRUE) +
      geom_contour(color = "black", alpha = 0.5, bins = 10) +
      coord_equal(expand = FALSE) +
      labs(
        x = "Downwind x (km)", y = "Crosswind y (km)",
        title = sprintf("PG plume at t = %.1f min (x_rel = x − U·t)", input$t_min),
        fill = if (input$pg_log) "log10(C) [g m⁻³]" else "C [g m⁻³]"
      ) +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(face = "bold"), legend.position = "right") +
      viridis::scale_fill_viridis(na.value = NA)
  })
  
  # Centerline slice at current t (y=0)
  output$pg_centerline <- renderPlot({
    xmax <- input$xmax_km * 1000
    dx   <- input$res_m
    xs <- seq(0, xmax, by = dx)
    t_sec <- input$t_min * 60
    x0 <- input$pg_U * t_sec
    x_rel <- xs - x0
    C <- plume_conc_pg(Q = qp_gps(), U = input$pg_U, H_eff = input$Hs,
                       x_rel = x_rel, y = rep(0, length(xs)), cls = input$pg_cls)
    df <- tibble(x = xs/1000, C = C)
    ggplot(df, aes(x, C)) +
      geom_line() +
      geom_vline(xintercept = input$focusX/1000, linetype = "dashed") +
      labs(x = "Downwind x (km)", y = "C at y=0 [g m⁻³]",
           title = "Centerline slice (temporal PG)") +
      theme_bw()
  })
  
  # Time series at chosen focusX (uses moving-frame arrival)
  output$pg_time_focus <- renderPlot({
    tgrid <- seq(input$t_window_min[1], input$t_window_min[2], by = input$dt_s/60)
    x_f <- input$focusX
    Cts <- map_dbl(tgrid, function(tmin) {
      x_rel <- x_f - input$pg_U * (tmin*60)
      plume_conc_pg(qp_gps(), input$pg_U, input$Hs, x_rel = x_rel, y = 0, cls = input$pg_cls)
    })
    df <- tibble(t = tgrid, C = Cts)
    ggplot(df, aes(t, C)) +
      geom_line() +
      geom_vline(xintercept = input$t_min, linetype = "dashed") +
      labs(x = "Time (min)", y = "C at focus X [g m⁻³]",
           title = sprintf("At focus X = %d m", input$focusX)) +
      theme_bw()
  })
}

shinyApp(ui, server)
