## Conceptual figures for paper 'Beyond Bergmann and Allen: A Unified Framework for Thermoregulatory Morphology'
# 

# Libraries ---------------------------------------------------------------
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(glue)
library(cowplot)
library(ggplot2)
library(png)
library(grid)
ggplot2::theme_set(theme_cowplot())

# Conceptual model -------------------------------------------------------
## Create conceptual model figure using DiagrammeR

# Function to paste pieces of figure together and depict using grViz(). I will create 'modules' that will be added to the base graph to build figure up with increasing complexity 
make_figure <- function(extra = "") {
  
  graph_txt <- glue("
  {base_graph}

  {extra}

  }}
  ")
  
  grViz(graph_txt)
}

# >Base graph -------------------------------------------------------------
# Base graph + interaction nodes
base_graph <- "
digraph conceptual_model {

graph [
  layout = neato,
  overlap = false,
  splines = true,
  outputorder = edgesfirst
]

node [
  fontname = Helvetica,
  fontsize = 18,
  color = black
]

edge [
  fontname = Helvetica,
  color = black,
  arrowsize = 1.2,
  penwidth = 1.8
]

# -------------------------
# Core nodes
# -------------------------

Temp [
  label = 'Temperature',
  shape = rectangle,
  width = 2.4,
  height = 0.9,
  pos = '0,1!'
]

Size [
  label = 'Size',
  shape = rectangle,
  width = 1.8,
  height = 0.8,
  pos = '5.8,2.3!'
]

Shape [
  label = 'Shape',
  shape = rectangle,
  width = 1.8,
  height = 0.8,
  pos = '5.8,-0.3!'
]

# -------------------------
# Core pathways
# -------------------------

Temp -> Size

Temp -> Shape

Size -> Shape [
  penwidth = 2.3
]

# -------------------------
# Interaction nodes
# -------------------------

Int1 [
  label = '',
  shape = point,
  width = 0.01,
  height = 0.01,
  style = invis,
  pos = '3.6, 0.22!'
]

Int2 [
  label = '',
  shape = point,
  width = 0.01,
  height = 0.01,
  style = invis,
  pos = '3.6, 1.8!'
]
"

# >Interactions -----------------------------------------------------------
# Interaction paths
interaction_paths <- "
# Interaction moderation

Size -> Int1 [
  style = dashed,
  color = gray40
]

Shape -> Int2 [
  style = dashed,
  color = gray40
]
"
interaction_grViz <- make_figure(interaction_paths)

# >Functional constraints -------------------------------------------------
# Functional constraints node + interaction paths
constraints_module <- "
# Functional constraints
FC [
  label = 'Functional\\nconstraints',
  shape = ellipse,
  width = 2,
  height = 1,
  pos = '2.5,1!'
]

FC -> Int1 [
  style = dashed,
  color = gray40
]

FC -> Int2 [
  style = dashed, 
  color = gray40
]
"
make_figure(constraints_module) 

# >Correlated responses ---------------------------------------------------
# Correlated responses between temperature and shape & temperature & size, as in complementarity 

# Want this to be curved
complementarity <- "
Int1 -> Int2 [
  dir = both,
  color = gray30,
  penwidth = 3,
  constraint = false,
  arrowsize = 0.9,
  minlen = 2,
  splines = curved
]
"

Baldwin_text <- glue(constraints_module) # , complementarity
Baldwin_grViz <- make_figure(Baldwin_text)
# Visualize
Baldwin_grViz

# >Full model -------------------------------------------------------------
# Paste together components of the model and export figure

Full_mod_text <- glue(interaction_paths, constraints_module, complementarity)
Full_mod_grViz <- make_figure(Full_mod_text)
# Visualize
Full_mod_grViz

# Export ------------------------------------------------------------------

svg <- export_svg(Full_mod_grViz)

export_graph <- function(g, file, width = 1200, height = 900){
  
  svg <- export_svg(g)
  
  rsvg_png(
    charToRaw(svg),
    file = file,
    width = width,
    height = height
  )
}

# Export as png
export_graph(interaction_grViz, "Figures/Conceptual_model_a.png")
export_graph(Baldwin_grViz, "Figures/Conceptual_model_b.png")

# Import and combine ------------------------------------------------------

img1 <- rasterGrob(readPNG("Figures/Conceptual_model_a.png"))
img2 <- rasterGrob(readPNG("Figures/Conceptual_model_b.png"))

p1 <- ggplot() +
  annotation_custom(img1, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void()

p2 <- ggplot() +
  annotation_custom(img2, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_void()

combined <- plot_grid(
  p1, p2,
  labels = c("a)", "b)"),
  ncol = 2,
  label_size = 14
)
# Visualized
combined 

# Save combined figure
ggsave(
  "Figures/Conceptual_model_panels.png",
  combined,
  width = 11.5,
  height = 4,
  dpi = 300,
  bg = "white"
)

stop()

# >Legend -----------------------------------------------------------------
## Working
Legend <- "
# -------------------------
# Legend
# -------------------------

subgraph cluster_legend {

  label = 'Legend'
  fontsize = 14
  fontname = Helvetica
  color = gray85
  style = rounded

  key1 [
    label = 'Direct effect',
    shape = plaintext,
    fontsize = 12
  ]

  key2 [
    label = 'Interaction / moderation',
    shape = plaintext,
    fontsize = 12
  ]

  solidA [
    label = '',
    shape = point,
    width = 0.01
  ]

  solidB [
    label = '',
    shape = point,
    width = 0.01
  ]

  dashA [
    label = '',
    shape = point,
    width = 0.01
  ]

  dashB [
    label = '',
    shape = point,
    width = 0.01
  ]

  solidA -> solidB [
    style = solid,
    penwidth = 2,
    arrowsize = 0.8
  ]

  dashA -> dashB [
    style = dashed,
    color = gray40,
    penwidth = 2,
    arrowsize = 0.8
  ]

  {rank = same; solidA; solidB; key1}
  {rank = same; dashA; dashB; key2}
}
"

Mod_legend <- paste0(Full_mod, Legend, "\n}")

make_figure(Mod_legend)

# Extras ------------------------------------------------------------------
# >IDAGs ------------------------------------------------------------------
library(dagitty)
library(ggdag)

## Shape as response
dag_shape <- dagify(
  size ~ temperature, 
  sizeXtemperature ~ size + temperature,
  shape ~ temperature + sizeXtemperature, # Ignore allometric effects of size for now
  exposure = c("sizeXtemperature"), 
  outcome = "shape",
  #latent = c() 
  labels = c("sizeXtemperature" = "size\nX\ntemperature")
)

tidy_dagitty(dag_shape, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none") +
  scale_fill_grey()

## Shape as response
dag_size_size <- dagify(
  size ~ temperature, 
  D_size_temperature ~ size, 
  exposure = c("size"), 
  outcome = "D_size_temperature"
  #latent = c() 
)

tidy_dagitty(dag_size_size, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")

## Temperature as response
dag_temp <- dagify(
  temperature ~ size + shape + sizeXshape, 
  sizeXshape ~ size + shape,
  shape ~ size,
  exposure = c("sizeXshape"), 
  outcome = "temperature"
  #latent = c() 
)

tidy_dagitty(dag_temp, layout = "fr") %>% 
  ggdag_status(text_col = "black",
               text = TRUE, 
               edge_type = "link_arc",
               node_size = 20,
               text_size = 3, 
               stylized = TRUE) + 
  theme_dag() + 
  guides(fill = "none", color = "none")

# Adjustment set
adjustmentSets(dag_shape, type = "minimal") 
