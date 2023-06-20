library(tidyverse)
library(magrittr)
library(liana)
library(ggplot2)
library(dplyr)

#please set the directory!
f_path <- "C:\\Users\\gocze\\OneDrive\\Biodata2\\NGS\\robustness data\\output\\liana_csvs"

rank_n <- 250
rank_list <- list('natmi.edge_specificity','connectome.weight_sc','logfc.logfc_comb','sca.LRscore')


results.df <- data.frame (Method=vector(),
                          CellsRemoved=vector(),
                          TPR=vector()) 

liana_df <- read.csv(paste(f_path,"\\liana_100.csv",sep=""))

#grand truth
gt_df <- data.frame()
env <- new.env()
col_name_list <- list('NATMI','Connectome','LogFC','SingleCellSignalR')


for (i in 1:length(rank_list)) {

expr <- parse(text = paste("gta_df <- liana_df[order(liana_df$", rank_list[[i]], ",decreasing = FALSE), ]", sep = ""))
eval(expr,env)

gta_df <- env$gta_df

gtm_df <- gta_df[1:rank_n, ]


gtm_df <- gtm_df %>%
  unite("LR_ID",
        c(source, target,ligand.complex, receptor.complex),
        remove = FALSE,
        sep = "_") %>%
  relocate("LR_ID", .after = last_col())

# Concatenate the columns horizontally
if (nrow(gt_df) == 0) {
  expr <- parse(text = paste("gt_df <- data.frame(", col_name_list[[i]], "=gtm_df[,c(\"LR_ID\")])", sep = ""))
  eval(expr,env)
  gt_df <- env$gt_df
} else {
  expr <- parse(text = paste("gt_df <- cbind(gt_df,data.frame(", col_name_list[[i]], "=gtm_df[,c(\"LR_ID\")]))", sep = ""))
  eval(expr,env)
  gt_df <- env$gt_df
}
}

perc <- list(0,20,40)

for (i in 1:length(col_name_list)) {
  for (j in 1:length(perc)) {
    if (perc[[j]] == 0){
      next_row_num = nrow(results.df) + 1
      print(next_row_num)
      results.df[next_row_num,1] <- col_name_list[[i]]
      results.df[next_row_num,2] <- perc[[j]]
      results.df[next_row_num,3] <- 1
    } else{
      for (k in 1:5) {
        sampled_df <- read.csv(paste(f_path,"\\uniform\\seu_",perc[[j]],"_u_",k,".csv",sep=""))
        s_df <- sampled_df[order(liana_df$natmi.edge_specificity,decreasing = FALSE), ]
        s_df <- s_df[1:rank_n, ]
        
        s_df <- s_df %>%
          unite("LR_ID",
                c(source, target,ligand.complex, receptor.complex),
                remove = FALSE,
                sep = "_") %>%
          relocate("LR_ID", .after = last_col())
        expr <- parse(text = paste("tp <- sum(s_df$LR_ID %in% gt_df$",col_name_list[[i]],")",sep=""))
        eval(expr,env)
        tp <- env$tp
        cp <- length(s_df$LR_ID)
        tpr <- tp / cp 
        next_row_num = nrow(results.df) + 1
        print(next_row_num)
        results.df[next_row_num,1] <- col_name_list[[i]]
        results.df[next_row_num,2] <- perc[[j]]
        results.df[next_row_num,3] <- tpr
      }
    }
  }
}

ylabel = "True Positive Rate"


results.df %>% 
  ggplot(aes(x = CellsRemoved,
             y = TPR,
             group = CellsRemoved,
             colour = Method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.4)+ 
  scale_y_continuous(breaks = seq(0, 1, 0.20), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 100, 20)) +
  # add text
  ylab(ylabel) +
  labs(subtitle = "Boxplot by Method.",
       color = "Method") +
  facet_grid(~Method, scales='free_x', space='free', switch="x") +
  theme_bw(base_size = 24) +
  theme(strip.text.x = element_text(angle = 90, face="bold", colour="white"),
        strip.background = element_rect(fill="darkgray"),
        legend.title = element_text(size = 28),
        legend.text = element_text(size = 25)
  ) +
  labs(colour=guide_legend(title="Method"),
       caption = NULL,
       title = NULL,
       subtitle = NULL) +
  xlab(str_glue("Cells Removed %")) +
  guides(shape = "none") +
  ggtitle("Subsampling with uniform")
