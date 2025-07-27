# ------------------------------------------------------------------
# admixtools2_interCladeAdmix_dxyDirectionTestt.R
#
# admixture graph estimation using admixtools2
# 2025/04/13, 
# Author: Shingo Fujimoto
# ------------------------------------------------------------------
library(admixtools)
library(ggplot2)
library(tidyverse)
library(igraph)

# Calculate and plot the admixturegraphs
# -------------------------------------------------------------------
qpGraphPlot <-  function(graphName)
{
  tmpScenario <- graphName
  graph_el <- parse_qpgraph_graphfile(tmpScenario, igraph = TRUE)
  fit <- qpgraph(F2_blocks, graph_el)
  graphOutput <-  plot_graph(fit$edges)
  graphOutput <- graphOutput + annotate("text", x=0, y=0, label = paste0("score = ", round(fit$score,2)))
  ggsave(file = paste0(graphName, ".pdf"), plot = graphOutput, dpi = 300, width = 6, height = 6)
  return(graphOutput)
}

qpGraphModel <-  function(graphName)
{
  tmpScenario <- graphName
  graph_el <- parse_qpgraph_graphfile(tmpScenario, igraph = TRUE)
  fit <- qpgraph(F2_blocks, graph_el, return_fstats = TRUE)
  return(fit)
}

# Admixturegraph in EJ (iso, mry), WJ, WK, Osak, HSOk 
# -------------------------------------------------------------------
snpDataPrefix <- "../genotype/EJsep_WJWK_Osak_HSOK/medakaWGS_HSOK_sak_lat_33samples"
popList <- c("eKor", "wKyu", "westJ", "sak", "iso", "mry")
outDIRname <- paste0('../genotype/admixtools2_medakaWGS_', paste(popList, collapse = '_'))

#F2blocks <- extract_f2(snpDataPrefix, outdir = outDIRname, pops = popList, n_cores = 4, auto_only = FALSE, overwrite = TRUE)
# import the pre calculated F2
F2_blocks <- f2_from_precomp(outDIRname)
#F2_blocks <- f2_from_precomp("../genotype/admixtools2_medakaWGS_eKor_wKyu_westJ_eastJ_sak_iso_mry") # import the pre calculated F2

# Estimate the admixturegraph in each scenario 
# noadmixture
qpGraphPlot("../qpgraph/admixGraph/medakaWGS_HSOK_sak_westK_westJ_eastJ__scenarioa")
# admixture, EJ(iso)->Osak
qpGraphPlot("../qpgraph/admixGraph/medakaWGS_HSOK_sak_westK_westJ_eastJ__scenariob")
# admixture, EJ->WJ(mry), WKyu->WJ
qpGraphPlot("../qpgraph/admixGraph/medakaWGS_HSOK_sak_westK_westJ_eastJ__scenarioc")
# admixture, Osak&WK -> WJ
qpGraphPlot("../qpgraph/admixGraph/medakaWGS_HSOK_sak_westK_westJ_eastJ__scenario5-d") # admixture, Osak&WK -> WJ
qpGraphPlot("../qpgraph/admixGraph/medakaWGS_HSOK_sak_westK_westJ_eastJ__scenario5-e")
qpGraphPlot("../qpgraph/admixGraph/medakaWGS_HSOK_sak_westK_westJ_eastJ__scenario5-f")
