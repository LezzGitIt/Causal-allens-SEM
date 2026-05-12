## Conceptual figures for paper 'Beyond Bergmann and Allen: A Unified Framework for Thermoregulatory Morphology'

# Libraries ---------------------------------------------------------------
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(glue)
ggplot2::theme_set(theme_cowplot())

# Conceptual model -------------------------------------------------------

# Function to paste new 'modules' to base graph and depict using grViz()
make_figure <- function(extra = "") {
  
  graph_txt <- glue("
  {base_graph}

  {extra}

  }}
  ")
  
  grViz(graph_txt)
}


# >Base graph -------------------------------------------------------------
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
make_figure(interaction_paths)


# >Functional constraints -------------------------------------------------
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

# >Full model -------------------------------------------------------------
# Temperature's influence on body size and shape is influenced by functional constraints, and by size and shape themselves. Size influences shape through allometry. Solid lines represent direct effects and dashed lines represent interactions.

Full_mod_text <- glue(interaction_paths, constraints_module)
Full_mod_grViz <- make_figure(Full_mod)
svg <- export_svg(Full_mod_grViz)

# Export as a png
rsvg_png(
  charToRaw(svg),
  file = "Figures/Conceptual_model.png",
  width = 850,
  height = 450
)

stop()

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
