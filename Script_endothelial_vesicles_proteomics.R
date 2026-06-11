
library(ggplot2)
library(BiocManager)
library(DESeq2)
library(dplyr)
library(edgeR)
library(sva)
library(tools)
library(dendextend)
library(limma)
library(preprocessCore)
library(rhdf5)
library(reshape)
library(EnhancedVolcano)
library(pheatmap)
library(devtools)
library(readr)


library(tidyverse)
df

df<-read_tsv('/Users/penri/OneDrive/Desktop/ProjetoUsp/AnaliseLeticia/Proteomicas-20240913/Fragpipe2/Data/report.pg_matrix.tsv')
df[df$Genes %in% "PLA2G6", ]

#TRANSFORMING TO DATAFRAME
df<-as.data.frame(df)


# Rename columns
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_2H.mzML"] <- "Heme2"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_1T.mzML"] <- "Tnf1"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_2C.mzML"] <- "Controle2"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_4T.mzML"] <- "Tnf4"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_4H.mzML"] <- "Heme4"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_5C.mzML"] <- "Controle5"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_3T.mzML"] <- "Tnf3"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_1C.mzML"] <- "Controle1"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_5T.mzML"] <- "Tnf5"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_5H.mzML"] <- "Heme5"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_4C.mzML"] <- "Controle4"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_1H.mzML"] <- "Heme1"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_3C.mzML"] <- "Controle3"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_2T.mzML"] <- "Tnf2"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_3H.mzML"] <- "Heme3"


df<-df[,-3]
df<-df[,-4]

df[is.na(df)]<-0


#Casting all columns to integer

df$Heme1<-as.integer(df$Heme1)
df$Heme2<-as.integer(df$Heme2)
df$Heme3<-as.integer(df$Heme3)
df$Heme4<-as.integer(df$Heme4)
df$Heme5<-as.integer(df$Heme5)
df$Controle1<-as.integer(df$Controle1)
df$Controle2<-as.integer(df$Controle2)
df$Controle3<-as.integer(df$Controle3)
df$Controle4<-as.integer(df$Controle4)
df$Controle5<-as.integer(df$Controle5)
df$Tnf1<-as.integer(df$Tnf1)
df$Tnf2<-as.integer(df$Tnf2)
df$Tnf3<-as.integer(df$Tnf3)
df$Tnf4<-as.integer(df$Tnf4)
df$Tnf5<-as.integer(df$Tnf5)

#Grouping by genes, and making then unique
df["Controle1"][is.na(df["Controle1"])] <- 0
df["Controle2"][is.na(df["Controle2"])] <- 0
df["Controle3"][is.na(df["Controle3"])] <- 0
df["Controle4"][is.na(df["Controle4"])] <- 0
df["Controle5"][is.na(df["Controle5"])] <- 0
df["Heme1"][is.na(df["Heme1"])] <- 0
df["Heme2"][is.na(df["Heme2"])] <- 0
df["Heme3"][is.na(df["Heme3"])] <- 0
df["Heme4"][is.na(df["Heme4"])] <- 0
df["Heme5"][is.na(df["Heme5"])] <- 0
df["Tnf1"][is.na(df["Tnf1"])] <- 0
df["Tnf2"][is.na(df["Tnf2"])] <- 0
df["Tnf3"][is.na(df["Tnf3"])] <- 0
df["Tnf4"][is.na(df["Tnf4"])] <- 0
df["Tnf5"][is.na(df["Tnf5"])] <- 0
df<-df[,-1]
df<-df[,-1]





#GROUPING BY
df<-as.data.frame(df)
library(dplyr)
df1<-group_by(df,Genes)

summarise(df, sum(Genes))

df1 <- df %>%
  group_by(Genes) %>%
  summarise(
    Controle1 = sum(Controle1),
    Controle2 = sum(Controle2),
    Controle3 = sum(Controle3),
    Controle4 = sum(Controle4),
    Controle5 = sum(Controle5),
    Heme1 = sum(Heme1),
    Heme2 = sum(Heme2),
    Heme3 = sum(Heme3),
    Heme4 = sum(Heme4),
    Heme5 = sum(Heme5),
    Tnf1 = sum(Tnf1),
    Tnf2 = sum(Tnf2),
    Tnf3 = sum(Tnf3),
    Tnf4 = sum(Tnf4),
    Tnf5 = sum(Tnf5)
  )


df1<- as.data.frame(df1)
rownames(df1)<- df1$Genes

df1<-df1[,-1]




#creating coldata for Deseq Analysis
coldata_dds<-as.data.frame(colnames(df1))
coldata_dds<-as.data.frame(coldata_dds)
coldata_dds$condition <-c('Controle','Controle','Controle',
                          'Controle','Controle' ,'Heme' ,
                          'Heme','Heme' ,'Heme',
                          'Heme','Tnf' ,'Tnf',
                          'Tnf','Tnf' ,'Tnf')


df1$Heme1<-as.integer(df1$Heme1)
df1$Heme2<-as.integer(df1$Heme2)
df1$Heme3<-as.integer(df1$Heme3)
df1$Heme4<-as.integer(df1$Heme4)
df1$Heme5<-as.integer(df1$Heme5)
df1$Controle1<-as.integer(df1$Controle1)
df1$Controle2<-as.integer(df1$Controle2)
df1$Controle3<-as.integer(df1$Controle3)
df1$Controle4<-as.integer(df1$Controle4)
df1$Controle5<-as.integer(df1$Controle5)
df1$Tnf1<-as.integer(df1$Tnf1)
df1$Tnf2<-as.integer(df1$Tnf2)
df1$Tnf3<-as.integer(df1$Tnf3)
df1$Tnf4<-as.integer(df1$Tnf4)
df1$Tnf5<-as.integer(df1$Tnf5)

df1



# DATA NNORMALIZATION
#BiocManager::install("vsn")
library(vsn)
## --------------------------------------------------------------------------
### VARIANCE STABILIZING NORMALIZATION
vsn_fit <- vsn2(as.matrix(df1))
df_vsn_norm <- predict(vsn_fit, newdata = as.matrix(df1))


#par(mfrow=c(1, 2))
#meanSdPlot(as.matrix(df1), main="Antes da VSN")

# Depois da normalização
#meanSdPlot(df_vsn_norm, main="Depois da VSN")
## --------------------------------------------------------------------------

#Quantile NORMALIZATION
library(preprocessCore)


# Normalização de quantis
df_quant_norm <- normalize.quantiles(as.matrix(df1))

# Converter para um data frame e adicionar nomes de genes
df_quant_norm <- as.data.frame(df_quant_norm)
rownames(df_quant_norm) <- rownames(df1)

colnames(df_quant_norm) <- colnames(df1)

#df1$Genes<-rownames(df1)
df1[df1$Genes %in% "GAPDH", ]

##------------------------------------------------------------

df_vsn_norm <- df_vsn_norm[apply(df_vsn_norm, 1, function(row) all(row >= -140)), ]

#df1 <- df1[apply(df1, 1, function(row) all(row <= ylim)), ]

set.seed(123)


#ylim<-100000000
boxplot(df_quant_norm,
        main = "Quantile normalization",
        ylim=c(0,ylim),
        xlab = "Samples", ylab = "Normalized Counts",
        col = "lightgray", frame = FALSE)


#df_quant_norm <- df_quant_norm[apply(df_quant_norm, 1, function(row) all(row <= ylim)), ]


df_melt_quant<-melt(df_quant_norm)
dados <- data.frame(valores = df_melt$value)

# Criar gráfico de densidade com personalização
p<-ggplot(dados, aes(x = valores)) +
  geom_density(fill = "lightblue", color = "blue", alpha = 0.7, size = 1.2) +  # Personalizar cores e espessura
  labs(title = "Gráfico de Densidade Personalizado", x = "Valores", y = "Densidade") +
  theme_minimal()  # Aplicar um tema limpo
p
df_melt_vsn<-melt(df1)
dados <- data.frame(valores = df_melt_vsn$value)

p<-ggplot(dados, aes(x = valores)) +
  geom_density(fill = "lightblue", color = "blue", alpha = 0.7, size = 1.2) +  # Personalizar cores e espessura
  labs(title = "Gráfico de Densidade Personalizado", x = "Valores", y = "Densidade") +
  theme_minimal()  # Aplicar um tema limpo
p




## --------------------------------------------------------------------------


df_quant_norm$protein_name <-rownames(df_quant_norm)

df_melt<-as.data.frame(melt(df_quant_norm))
df_melt$group<-df_melt$variable 


df_melt$group <- gsub("Controle1", "Controle", df_melt$group)
df_melt$group <- gsub("Controle2", "Controle", df_melt$group)
df_melt$group <- gsub("Controle3", "Controle", df_melt$group)
df_melt$group <- gsub("Controle4", "Controle", df_melt$group)
df_melt$group <- gsub("Controle5", "Controle", df_melt$group)
df_melt$group <- gsub("Heme1", "Heme", df_melt$group)
df_melt$group <- gsub("Heme2", "Heme", df_melt$group)
df_melt$group <- gsub("Heme3", "Heme", df_melt$group)
df_melt$group <- gsub("Heme4", "Heme", df_melt$group)
df_melt$group <- gsub("Heme5", "Heme", df_melt$group)
df_melt$group <- gsub("Tnf1", "Tnf", df_melt$group)
df_melt$group <- gsub("Tnf2", "Tnf", df_melt$group)
df_melt$group <- gsub("Tnf3", "Tnf", df_melt$group)
df_melt$group <- gsub("Tnf4", "Tnf", df_melt$group)
df_melt$group <- gsub("Tnf5", "Tnf", df_melt$group)





df_norm<-as.data.frame(df_quant_norm)

df_norm$Heme1<-as.integer(df_norm$Heme1)
df_norm$Heme2<-as.integer(df_norm$Heme2)
df_norm$Heme3<-as.integer(df_norm$Heme3)
df_norm$Heme4<-as.integer(df_norm$Heme4)
df_norm$Heme5<-as.integer(df_norm$Heme5)
df_norm$Controle1<-as.integer(df_norm$Controle1)
df_norm$Controle2<-as.integer(df_norm$Controle2)
df_norm$Controle3<-as.integer(df_norm$Controle3)
df_norm$Controle4<-as.integer(df_norm$Controle4)
df_norm$Controle5<-as.integer(df_norm$Controle5)
df_norm$Tnf1<-as.integer(df_norm$Tnf1)
df_norm$Tnf2<-as.integer(df_norm$Tnf2)
df_norm$Tnf3<-as.integer(df_norm$Tnf3)
df_norm$Tnf4<-as.integer(df_norm$Tnf4)
df_norm$Tnf5<-as.integer(df_norm$Tnf5)

###-------------------------------------------------------------------------------------------------------------------------------
#Grafico de boxplot para algumas proteinas específicas


# Vamos supor que seus dados estão no data frame 'df'
library(ggplot2)

# Ajustar os dados para o formato correto
df_melt$protein_name <- as.factor(df_melt$protein_name)  # Certifique-se de que o nome da proteína é um fator
df_melt$group <- as.factor(df_melt$group)  # Garantir que o grupo também seja um fator

df_melt <- df_melt %>%
  select(-variable)

df_genes<-df_melt
df_melt[df_melt$protein_name %in% "GAPDH", ]


df_genes$value<- df_genes$value+1


# Exemplo de uma lista de genes de interesse
genes_de_interesse1 <- c("PRCP", "SRPX","EPHB2","COL15A1","PGF","FLT1","CCN1")  # Adicione os genes de interesse

genes_de_interesse2 <- c("RTN4", "RPS7","RPL10","FAP","","","")  # Adicione os genes de interesse

genes_de_interesse3 <- c("DKK1","TXNDC12","ANXA1","","","","")  # Adicione os genes de interesse

genes_de_interesse4 <- c("S100A8", "S100A9","GAPDH","PLA2G6","ASAH1","PPARG","GRN")  # Adicione os genes de interesse

# Filtrando os dados para manter apenas os genes da lista
df_genes <- df_genes[df_genes$protein_name %in% genes_de_interesse1, ]
df_genes$value2<-log2(df_genes$value)

## Teste normalidade

shapiro.test(df_genes$value2) ### genes de interesse 01 segue distribuicao normal

View(df_genes)




# Teste de Tukey
tukey1 <- TukeyHSD(anova1)
print(tukey1)

# Extraia a tabela de comparações do teste de Tukey para o fator "Sample"
df_tukey_resultBCA <- as.data.frame(tukey_resultBCA$Sample)

# Opcional: se quiser adicionar uma coluna com os pares de comparação (como "Heme-Control")
df_tukey_resultBCA$comparison <- rownames(df_tukey_resultBCA)

# Exiba o data frame
print(df_tukey_resultBCA)

# Extrair uma linha específica, como a comparação "Heme-Control"
Heme_ControlBCA <- df_tukey_resultBCA["HEME-CTRL", ]

TNF_ControlBCA <- df_tukey_resultBCA["TNF-CTRL", ]

TNF_HemeBCA <- df_tukey_resultBCA["TNF-HEME", ]



# install.packages("gridExtra")

df_summary <- df_genes %>%
  group_by(protein_name, group) %>%
  summarise(mean_value = mean(value2, na.rm = TRUE),
            sd_value = sd(value2, na.rm = TRUE)) %>%
  ungroup()





ggplot() +
  # Pontos com preenchimento e contorno diferenciados
  geom_point(
    data = df_genes,
    aes(
      x = protein_name,
      y = value2,
      fill = group,  # Preenchimento de acordo com o grupo
      color = group  # Contorno de acordo com o grupo
    ),
    size = 4,
    position = position_dodge(width = 0.8),
    shape = 21,  # Shape que permite preenchimento e contorno
    stroke = 0.8  # Espessura do contorno
  ) +
  
  # Barras de erro com cores dos grupos
  geom_errorbar(
    data = df_summary,
    aes(
      x = protein_name,
      ymin = mean_value - sd_value,
      ymax = mean_value + sd_value,
      color = group  # Cor das barras de erro
    ),
    width = 0.1,
    size = 0.4,
    position = position_dodge(width = 0.8)
  ) +
  
  # Rótulos e títulos
  labs(
    title = "Protein Abundance with Error Bars",
    x = "Proteins",
    y = "Normalized Protein Abundance (log2)",
    color = "",
    fill = ""
  ) +
  
  # Escalas de cor
  scale_fill_manual(values = c("#F4A6C8", "#A0D8B3", "#7BA9E3")) +  # Preenchimento
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +  # Contorno
  
  # Tema e personalização
  theme_minimal() +
  ylim(0, 50) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotaciona os nomes das proteínas no eixo X
    axis.title.x = element_text(face = "bold", size = 14),  # Título do eixo X
    axis.title.y = element_text(face = "bold", size = 14),  # Título do eixo Y
    axis.text.y = element_text(size = 12),  # Tamanho dos textos no eixo Y
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black", size = 1.2)
  )








###--------------------------------------------------------------------------------------------------
# Ajustar os dados para o formato correto
df_melt$protein_name <- as.factor(df_melt$protein_name)  # Certifique-se de que o nome da proteína é um fator
df_melt$group <- as.factor(df_melt$group)  # Garantir que o grupo também seja um fator

df_melt <- df_melt %>%
  select(-variable)

df_genes<-df_melt
df_melt[df_melt$protein_name %in% "GAPDH", ]


df_genes$value<- df_genes$value+1


# Exemplo de uma lista de genes de interesse
genes_de_interesse1 <- c("PRCP", "SRPX","EPHB2","COL15A1","PGF","FLT1","CCN1")  # Adicione os genes de interesse

genes_de_interesse2 <- c("RTN4", "RPS7","RPL10","FAP","","","")  # Adicione os genes de interesse

genes_de_interesse3 <- c("DKK1","TXNDC12","ANXA1","","","","")  # Adicione os genes de interesse

genes_de_interesse4 <- c("S100A8", "S100A9","GAPDH","PLA2G6","ASAH1","PPARG","GRN")  # Adicione os genes de interesse

# Filtrando os dados para manter apenas os genes da lista
df_genes <- df_genes[df_genes$protein_name %in% genes_de_interesse_ferroptosis, ]
df_genes$value2<-log2(df_genes$value)


# install.packages("gridExtra")

df_summary <- df_genes %>%
  group_by(protein_name, group) %>%
  summarise(mean_value = mean(value2, na.rm = TRUE),
            sd_value = sd(value2, na.rm = TRUE)) %>%
  ungroup()





ggplot() +
  # Pontos com preenchimento e contorno diferenciados
  geom_point(
    data = df_genes,
    aes(
      x = protein_name,
      y = value2,
      fill = group,  # Preenchimento de acordo com o grupo
      color = group  # Contorno de acordo com o grupo
    ),
    size = 4,
    position = position_dodge(width = 0.8),
    shape = 21,  # Shape que permite preenchimento e contorno
    stroke = 0.8  # Espessura do contorno
  ) +
  
  # Barras de erro com cores dos grupos
  geom_errorbar(
    data = df_summary,
    aes(
      x = protein_name,
      ymin = mean_value - sd_value,
      ymax = mean_value + sd_value,
      color = group  # Cor das barras de erro
    ),
    width = 0.1,
    size = 0.4,
    position = position_dodge(width = 0.8)
  ) +
  
  # Rótulos e títulos
  labs(
    title = "Protein Abundance with Error Bars",
    x = "Proteins",
    y = "Normalized Protein Abundance (log2)",
    color = "",
    fill = ""
  ) +
  
  # Escalas de cor
  scale_fill_manual(values = c("#F4A6C8", "#A0D8B3", "#7BA9E3")) +  # Preenchimento
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +  # Contorno
  
  # Tema e personalização
  theme_minimal() +
  ylim(0, 50) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotaciona os nomes das proteínas no eixo X
    axis.title.x = element_text(face = "bold", size = 14),  # Título do eixo X
    axis.title.y = element_text(face = "bold", size = 14),  # Título do eixo Y
    axis.text.y = element_text(size = 12),  # Tamanho dos textos no eixo Y
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black", size = 1.2)
  )




























###--------------------------------------------------------------------------------------------------
# Ajustar os dados para o formato correto
df_melt$protein_name <- as.factor(df_melt$protein_name)  # Certifique-se de que o nome da proteína é um fator
df_melt$group <- as.factor(df_melt$group)  # Garantir que o grupo também seja um fator

df_melt <- df_melt %>%
  select(-variable)

df_genes<-df_melt
df_melt[df_melt$protein_name %in% "GAPDH", ]


df_genes$value<- df_genes$value+1


# Exemplo de uma lista de genes de interesse
genes_de_interesse1 <- c("PRCP", "SRPX","EPHB2","COL15A1","PGF","FLT1","CCN1")  # Adicione os genes de interesse

genes_de_interesse2 <- c("RTN4", "RPS7","RPL10","FAP","","","")  # Adicione os genes de interesse

genes_de_interesse3 <- c("DKK1","TXNDC12","ANXA1","","","","")  # Adicione os genes de interesse

genes_de_interesse4 <- c("S100A8", "S100A9","GAPDH","PLA2G6","ASAH1","PPARG","GRN")  # Adicione os genes de interesse

# Filtrando os dados para manter apenas os genes da lista
df_genes <- df_genes[df_genes$protein_name %in% genes_de_interesse_ferroptosis, ]
df_genes$value2<-log2(df_genes$value)


# install.packages("gridExtra")

df_summary <- df_genes %>%
  group_by(protein_name, group) %>%
  summarise(mean_value = mean(value2, na.rm = TRUE),
            sd_value = sd(value2, na.rm = TRUE)) %>%
  ungroup()





ggplot() +
  # Pontos com preenchimento e contorno diferenciados
  geom_point(
    data = df_genes,
    aes(
      x = protein_name,
      y = value2,
      fill = group,  # Preenchimento de acordo com o grupo
      color = group  # Contorno de acordo com o grupo
    ),
    size = 4,
    position = position_dodge(width = 0.8),
    shape = 21,  # Shape que permite preenchimento e contorno
    stroke = 0.8  # Espessura do contorno
  ) +
  
  # Barras de erro com cores dos grupos
  geom_errorbar(
    data = df_summary,
    aes(
      x = protein_name,
      ymin = mean_value - sd_value,
      ymax = mean_value + sd_value,
      color = group  # Cor das barras de erro
    ),
    width = 0.1,
    size = 0.4,
    position = position_dodge(width = 0.8)
  ) +
  
  # Rótulos e títulos
  labs(
    title = "Protein Abundance with Error Bars",
    x = "Proteins",
    y = "Normalized Protein Abundance (log2)",
    color = "",
    fill = ""
  ) +
  
  # Escalas de cor
  scale_fill_manual(values = c("#F4A6C8", "#A0D8B3", "#7BA9E3")) +  # Preenchimento
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +  # Contorno
  
  # Tema e personalização
  theme_minimal() +
  ylim(0, 50) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotaciona os nomes das proteínas no eixo X
    axis.title.x = element_text(face = "bold", size = 14),  # Título do eixo X
    axis.title.y = element_text(face = "bold", size = 14),  # Título do eixo Y
    axis.text.y = element_text(size = 12),  # Tamanho dos textos no eixo Y
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black", size = 1.2)
  )
























###--------------------------------------------------------------------------------------------------
# Ajustar os dados para o formato correto
df_melt$protein_name <- as.factor(df_melt$protein_name)  # Certifique-se de que o nome da proteína é um fator
df_melt$group <- as.factor(df_melt$group)  # Garantir que o grupo também seja um fator

df_melt <- df_melt %>%
  select(-variable)

df_genes<-df_melt
df_melt[df_melt$protein_name %in% "GAPDH", ]


df_genes$value<- df_genes$value+1


# Exemplo de uma lista de genes de interesse
genes_de_interesse1 <- c("PRCP", "SRPX","EPHB2","COL15A1","PGF","FLT1","CCN1")  # Adicione os genes de interesse

genes_de_interesse2 <- c("RTN4", "RPS7","RPL10","FAP","","","")  # Adicione os genes de interesse

genes_de_interesse3 <- c("DKK1","TXNDC12","ANXA1","","","","")  # Adicione os genes de interesse

genes_de_interesse4 <- c("S100A8", "S100A9","GAPDH","PLA2G6","ASAH1","PPARG","GRN")  # Adicione os genes de interesse

# Filtrando os dados para manter apenas os genes da lista
df_genes <- df_genes[df_genes$protein_name %in% genes_de_interesse_ferroptosis, ]
df_genes$value2<-log2(df_genes$value)


# install.packages("gridExtra")

df_summary <- df_genes %>%
  group_by(protein_name, group) %>%
  summarise(mean_value = mean(value2, na.rm = TRUE),
            sd_value = sd(value2, na.rm = TRUE)) %>%
  ungroup()





ggplot() +
  # Pontos com preenchimento e contorno diferenciados
  geom_point(
    data = df_genes,
    aes(
      x = protein_name,
      y = value2,
      fill = group,  # Preenchimento de acordo com o grupo
      color = group  # Contorno de acordo com o grupo
    ),
    size = 4,
    position = position_dodge(width = 0.8),
    shape = 21,  # Shape que permite preenchimento e contorno
    stroke = 0.8  # Espessura do contorno
  ) +
  
  # Barras de erro com cores dos grupos
  geom_errorbar(
    data = df_summary,
    aes(
      x = protein_name,
      ymin = mean_value - sd_value,
      ymax = mean_value + sd_value,
      color = group  # Cor das barras de erro
    ),
    width = 0.1,
    size = 0.4,
    position = position_dodge(width = 0.8)
  ) +
  
  # Rótulos e títulos
  labs(
    title = "Protein Abundance with Error Bars",
    x = "Proteins",
    y = "Normalized Protein Abundance (log2)",
    color = "",
    fill = ""
  ) +
  
  # Escalas de cor
  scale_fill_manual(values = c("#F4A6C8", "#A0D8B3", "#7BA9E3")) +  # Preenchimento
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +  # Contorno
  
  # Tema e personalização
  theme_minimal() +
  ylim(0, 50) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotaciona os nomes das proteínas no eixo X
    axis.title.x = element_text(face = "bold", size = 14),  # Título do eixo X
    axis.title.y = element_text(face = "bold", size = 14),  # Título do eixo Y
    axis.text.y = element_text(size = 12),  # Tamanho dos textos no eixo Y
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black", size = 1.2)
  )




















###-------------------------------------------------------------------------------------------------------------------------------


packageVersion("DESeq2")

#DESEQ2
dds<-DESeqDataSetFromMatrix(countData = as.matrix(df_norm),
                            colData = coldata_dds, 
                            design = ~condition)

dds<-DESeq(dds)

res<-results(dds, contrast=c('condition', 'Tnf','Controle'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)


ggplot(res, aes(x = log2FoldChange, y = negLogp)) +
  theme_bw() +  
  theme(
    plot.title = element_text(size = 35, face = "bold", hjust = 0.5),  # Centraliza o título
    axis.title = element_text(size = 45),  # Aumenta o tamanho do título dos eixos
    axis.text = element_text(size = 45),   # Aumenta o tamanho dos valores dos eixos
    panel.border = element_rect(size = 2.5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "none"
  ) +
  geom_point(aes(color = ifelse(negLogp > -log10(pvalue_threshold) & abs(log2FoldChange) > log2FC_threshold, "significant", "non-significant")), 
             alpha = 0.6) +
  scale_color_manual(values = c("significant" = "red", "non-significant" = "black")) +  
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = log2FC_threshold, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = -log2FC_threshold, linetype = "dashed", color = "blue") +
  geom_text(data = subset(res), 
            aes(label = genes), 
            vjust = 1.5, 
            size = 12, 
            check_overlap = TRUE) +
  labs(color = "") +
  xlim(-6,6)+
  ylim(0,5)+
  ggtitle("Control X TNF")

# ggsave('C:/Users/penri/OneDrive/Desktop/ProjetoUsp/AnaliseLeticia/Proteomicas-20240913/Fragpipe2/Graphs/Vulcano_plot/volcano_plot_Control_TNF_20241228.svg')
# library('svglite')





#png("C:/Users/penri/Downloads/Volcano_plot_Control_Heme_20241226.png", width = 1200, height = 800)  # Defina o caminho onde salvar
#draw(vp)
#dev.off()


#df_norm[rownames(df_norm) == "ACTB", ]



##----------------------------------------------------
# Install packages if not already installed
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}

if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}

# Load libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(grid)   # needed for gpar() and unit()

res<-results(dds, contrast=c('condition', 'Heme','Controle'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)

filtered_res_heme <- res[which((res$pvalue<=0.05) & abs(res$log2FoldChange) > 1), ]
#filtered_res_heme <- res[which(abs(res$log2FoldChange) > 1), ]

res<-results(dds, contrast=c('condition', 'Tnf','Controle'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)

filtered_res_tnf <- res[which((res$pvalue<=0.05) & abs(res$log2FoldChange) > 1), ]
#filtered_res_tnf <- res[which(abs(res$log2FoldChange) > 1), ]



res<-results(dds, contrast=c('condition', 'Tnf','Heme'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)

filtered_res_Heme_Tnf <- res[which((res$pvalue<=0.05) & abs(res$log2FoldChange) > 1), ]
#filtered_res_tnf <- res[which(abs(res$log2FoldChange) > 1), ]

dim(filtered_res_Heme_Tnf)


### AGORA VAMOS CRIAR UM MERGE PARA QUE APARECAM AS INFORMAOES DE AMBOS OS FILTROS
# Transformar os resultados em data frames, se necessário
filtered_res_tnf_df <- as.data.frame(filtered_res_tnf)
filtered_res_heme_df <- as.data.frame(filtered_res_heme)
filtered_res_heme_tnf_df <- as.data.frame(filtered_res_Heme_Tnf)

# Adicionar uma coluna identificadora para cada conjunto de resultados
filtered_res_tnf_df$condition <- "TNF"
filtered_res_heme_df$condition <- "Heme"
filtered_res_heme_tnf_df$condition<- "Heme_Tnf"

# Union das duas tabelas usando rbind()
merged_results <- rbind(filtered_res_tnf_df, filtered_res_heme_df,filtered_res_heme_tnf_df)

top_genes_tnf<-rownames(filtered_res_tnf_df)
top_genes_heme<-rownames(filtered_res_heme_df)
top_genes_heme_tnf<-rownames(filtered_res_heme_tnf_df)


top_genes<-as.data.frame(unique(merged_results$genes))
dim(top_genes)

#write.csv(top_genes , file = "top_genes.csv", row.names = FALSE)

library(dplyr)


df_norm$genes<-rownames(df_norm)

df_norm[df_norm$genes == 'FAP',]
df_norm<-df_norm[,1:15]

#norm_counts<-counts(dds)
norm_z<-df_norm[rownames(df_norm) %in% top_genes[,1], ]
norm_z<-norm_z[,1:15]

norm_z<-t(apply(norm_z,1,scale))
#colnames(norm_z)<-coldata_dds$`colnames(df1)`


#FAZER A DIFERENÇA ENTRE A MÉDIA ENTRE CONTROLE E OS TRATAMENTOS

#df_norm$CONTROL<-rowMeans(df_norm[, 1:5])

#df_norm<- df_norm/ df_norm$CONTROL
#df_norm$TNF <- rowMeans(df_norm[, 11:15])/ df_norm$CONTROL
#df_heatmap<-df_norm[,0:15]
#df_heatmap<-df_heatmap[rownames(df_heatmap) %in% top_genes[,1], ]
#df_heatmap<-t(apply(df_heatmap,1,scale))


ncol(df_heatmap)


df_heatmap<-norm_z
colnames(df_heatmap) <- c("CONTROL","CONTROL","CONTROL","CONTROL","CONTROL",
                          "HEME", "HEME", "HEME","HEME","HEME",
                          "TNF", "TNF", "TNF","TNF","TNF")


ht <- Heatmap(
  df_heatmap,
  name = "Z-score",                 # Legend title
  cluster_rows = TRUE,              # Cluster genes
  cluster_columns = FALSE,           # Keep sample order fixed
  column_labels = colnames(df_heatmap),
  column_names_gp = gpar(fontsize = 14),
  row_names_gp = gpar(fontsize = 10),
  col = colorRampPalette(brewer.pal(9, "Oranges"))(100),
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = NULL,
    legend_width = unit(6, "cm"),
    legend_height = unit(0.6, "cm"),
    direction = "horizontal",
    labels_gp = gpar(fontsize = 12),
    border = "black"
  ),
  column_title = "Differentially Expressed Genes",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_names_rot = 90              # Rotate column labels
)

ht
# png("C:/Users/penri/Downloads/heatmap_top_genes_20241226.png", width = 800, height = 1500)  # Defina o caminho onde salvar
draw(ht, heatmap_legend_side = "bottom")
dev.off()

# Obter o caminho do diretório "Downloads" para Windows
# downloads_path <- file.path(Sys.getenv("USERPROFILE"), "Downloads")
# downloads_path


#MONTAGEM DO HEATMAP
library(ComplexHeatmap)
packageVersion("ggplot2")

colnames(norm_z) <- c("CONTROL", "CONTROL", "CONTROL","CONTROL","CONTROL", 
                      "HEME", "HEME", "HEME","HEME","HEME",
                      "TNF", "TNF", "TNF","TNF","TNF")

library(RColorBrewer)

library(ComplexHeatmap)
library(circlize)  # for colorRamp2()
library(grid)
  

Heatmap(
  norm_z,
  name = "Z-score",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_labels = colnames(norm_z),
  column_names_gp = gpar(fontsize = 16),
  row_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    title = "Z-score",
    direction = "horizontal"
  )
)
dev.off()



res_heme<-as.data.frame(res_heme[rownames(res_heme) %in% (top_genes[,1]),])
res_tnf<-as.data.frame(res_tnf[rownames(res_tnf) %in% (top_genes[,1]),])



merged_results<-merge(res_heme, res_tnf,by="genes", all= T)

rownames(merged_results)<- merged_results$genes

temp<-merged_results[, c("log2FoldChange.x", "log2FoldChange.y")]
colnames(temp) <- c("Heme", "TNF")

Heatmap(temp,
        cluster_rows = T,
        cluster_columns = F,
        column_labels = colnames(temp),
        column_names_gp = gpar(fontsize = 16),
        row_names_gp = gpar(fontsize = 13))




###VENN DIAGRAMM
library(VennDiagram)



df<-read_tsv('/Users/penri/OneDrive/Desktop/ProjetoUsp/AnaliseLeticia/Proteomicas-20240913/Fragpipe2/Data/report.pg_matrix.tsv')



#TRANSFORMING TO DATAFRAME
df<-as.data.frame(df)


# Rename columns
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_2H.mzML"] <- "Heme2"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_1T.mzML"] <- "Tnf1"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_2C.mzML"] <- "Controle2"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_4T.mzML"] <- "Tnf4"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_4H.mzML"] <- "Heme4"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_5C.mzML"] <- "Controle5"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_3T.mzML"] <- "Tnf3"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_1C.mzML"] <- "Controle1"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_5T.mzML"] <- "Tnf5"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_5H.mzML"] <- "Heme5"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_4C.mzML"] <- "Controle4"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_1H.mzML"] <- "Heme1"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_3C.mzML"] <- "Controle3"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_2T.mzML"] <- "Tnf2"
names(df)[names(df) == "/project/yutaka/transcriptomeMetaAnalysis/rawData/leticiaFernandes/arquivosMZML/Leticia_08302024_3H.mzML"] <- "Heme3"


df<-df[,-3]
df<-df[,-4]

df[is.na(df)]<-0


#Casting all columns to integer

df$Heme1<-as.integer(df$Heme1)
df$Heme2<-as.integer(df$Heme2)
df$Heme3<-as.integer(df$Heme3)
df$Heme4<-as.integer(df$Heme4)
df$Heme5<-as.integer(df$Heme5)
df$Controle1<-as.integer(df$Controle1)
df$Controle2<-as.integer(df$Controle2)
df$Controle3<-as.integer(df$Controle3)
df$Controle4<-as.integer(df$Controle4)
df$Controle5<-as.integer(df$Controle5)
df$Tnf1<-as.integer(df$Tnf1)
df$Tnf2<-as.integer(df$Tnf2)
df$Tnf3<-as.integer(df$Tnf3)
df$Tnf4<-as.integer(df$Tnf4)
df$Tnf5<-as.integer(df$Tnf5)

#Grouping by genes, and making then unique
df["Controle1"][is.na(df["Controle1"])] <- 0
df["Controle2"][is.na(df["Controle2"])] <- 0
df["Controle3"][is.na(df["Controle3"])] <- 0
df["Controle4"][is.na(df["Controle4"])] <- 0
df["Controle5"][is.na(df["Controle5"])] <- 0
df["Heme1"][is.na(df["Heme1"])] <- 0
df["Heme2"][is.na(df["Heme2"])] <- 0
df["Heme3"][is.na(df["Heme3"])] <- 0
df["Heme4"][is.na(df["Heme4"])] <- 0
df["Heme5"][is.na(df["Heme5"])] <- 0
df["Tnf1"][is.na(df["Tnf1"])] <- 0
df["Tnf2"][is.na(df["Tnf2"])] <- 0
df["Tnf3"][is.na(df["Tnf3"])] <- 0
df["Tnf4"][is.na(df["Tnf4"])] <- 0
df["Tnf5"][is.na(df["Tnf5"])] <- 0
df<-df[,-1]
df<-df[,-1]




#GROUPING BY
df<-as.data.frame(df)
library(dplyr)
df1<-group_by(df,Genes)

summarise(df, sum(Genes))


df1 <- df %>%
  group_by(Genes) %>%
  summarise(
    Controle1 = sum(Controle1),
    Controle2 = sum(Controle2),
    Controle3 = sum(Controle3),
    Controle4 = sum(Controle4),
    Controle5 = sum(Controle5),
    Heme1 = sum(Heme1),
    Heme2 = sum(Heme2),
    Heme3 = sum(Heme3),
    Heme4 = sum(Heme4),
    Heme5 = sum(Heme5),
    Tnf1 = sum(Tnf1),
    Tnf2 = sum(Tnf2),
    Tnf3 = sum(Tnf3),
    Tnf4 = sum(Tnf4),
    Tnf5 = sum(Tnf5)
  )


df1<- as.data.frame(df1)
rownames(df1)<- df1$Genes

df1<-df1[,-1]




df1$Heme1<-as.integer(df1$Heme1)
df1$Heme2<-as.integer(df1$Heme2)
df1$Heme3<-as.integer(df1$Heme3)
df1$Heme4<-as.integer(df1$Heme4)
df1$Heme5<-as.integer(df1$Heme5)
df1$Controle1<-as.integer(df1$Controle1)
df1$Controle2<-as.integer(df1$Controle2)
df1$Controle3<-as.integer(df1$Controle3)
df1$Controle4<-as.integer(df1$Controle4)
df1$Controle5<-as.integer(df1$Controle5)
df1$Tnf1<-as.integer(df1$Tnf1)
df1$Tnf2<-as.integer(df1$Tnf2)
df1$Tnf3<-as.integer(df1$Tnf3)
df1$Tnf4<-as.integer(df1$Tnf4)
df1$Tnf5<-as.integer(df1$Tnf5)

dim(df1)


Control_sum<- as.data.frame(rowSums(df1[, c("Controle1", "Controle2", "Controle3", "Controle4", "Controle5")]))
Heme_sum <- as.data.frame(rowSums(df1[, c("Heme1", "Heme2", "Heme3", "Heme4", "Heme5")]))
Tnf_sum<- as.data.frame(rowSums(df1[, c("Tnf1", "Tnf2", "Tnf3", "Tnf4", "Tnf5")]))


Control_sum<-as.data.frame(rownames(Control_sum)[Control_sum>0])
Heme_sum<-as.data.frame(rownames(Heme_sum)[Heme_sum>0])
Tnf_sum<-as.data.frame(rownames(Tnf_sum)[Tnf_sum>0])


dim(Control_sum)
dim(Heme_sum)
dim(Tnf_sum)

options(max.print = 1600) 

proteins_control <- as.data.frame(unique(Control_sum[,1]))
proteins_heme <- as.data.frame(unique(Heme_sum[,1]))
proteins_tnf <- as.data.frame(unique(Tnf_sum[,1]))





# Criar o diagrama de Venn
venn.plot <- venn.diagram(
  x = list(
    Control = proteins_control,
    Heme = proteins_heme,
    TNF = proteins_tnf
  ),
  category.names = c("Control", "Heme", "TNF"),
  filename = NULL,
  output = TRUE,
  fill = c("blue", "red", "green"),
  alpha = 0.5,
  label.col = "white",
  cex = 1.5,
  cat.cex = 1.5,
  lwd = 2
)
grid.draw(venn.plot)

# Identificando proteínas exclusivas
unique_control <- setdiff(proteins_control, union(proteins_heme, proteins_tnf))
unique_heme <- setdiff(proteins_heme, union(proteins_control, proteins_tnf))
unique_tnf <- setdiff(proteins_tnf, union(proteins_control, proteins_heme))

# Exibir resultados
cat("Proteínas exclusivas do controle:\n", unique_control, "\n")
cat("Proteínas exclusivas do heme:\n", unique_heme, "\n")
cat("Proteínas exclusivas do TNF:\n", unique_tnf, "\n")














#GENE ANOTATION
res<-results(dds, contrast=c('condition', 'Heme','Controle'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)

filtered_res_heme <- res[which(res$pvalue < 0.05), ]
#filtered_res_heme <- res[which(abs(res$log2FoldChange) > 1), ]


res<-results(dds, contrast=c('condition', 'Tnf','Controle'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)

filtered_res_tnf <- res[which(res$pvalue < 0.05), ]




res<-results(dds, contrast=c('condition', 'Tnf','Heme'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)

filtered_res_tnf_Heme <- res[which(res$pvalue < 0.05), ]



a<-(genes = as.data.frame(rownames(filtered_res_heme)))
b<-(genes=as.data.frame(rownames(filtered_res_tnf)))
c<-(genes=as.data.frame(rownames(filtered_res_tnf_Heme)))

rbind(a,b)


#write.csv(df, "\\Users\\penri\\OneDrive\\Desktop\\ProjetoUsp\\AnaliseLeticia\\Proteomicas-20240913\\Fragpipe2\\Resultados_Path\\Análise_GO\\df.csv")







#######---------------------------------------------------------------------------
#Fazer o gráfico de Gene Ontology




#Genes que estão ligados com a coagulação sanguinea
top_genes<-data.frame(top_genes = c('GRN','FAP','LOX','GFER','LAMTOR5','PSME3','HS1BP3','CCAR2','NQO1','RTN4','TJP1',
                                    'PAK3','S100A8','S100A9','TXNDC12','DKK1','ANXA1','EDN1','RPS6','ASAH1','RPL10',
                                    'BCL2L13','RPS7','RPL11','FGA'))




dim(top_genes)[1]
dim(df1)
# Normalização de quantis
df_quant_norm <- normalize.quantiles(as.matrix(df1))

# Converter para um data frame e adicionar nomes de genes
df_quant_norm <- as.data.frame(df_quant_norm)
rownames(df_quant_norm) <- rownames(df1)

colnames(df_quant_norm) <- colnames(df1)


df_norm <- df_quant_norm
dim(df_norm)



# Calcular a média de expressão para cada gene entre as amostras
df_norm$CONTROL <- rowMeans(df_norm[, 1:5])
df_norm$HEME <- rowMeans(df_norm[, 6:10])
df_norm$TNF <- rowMeans(df_norm[, 11:15])


df_norm$gene<-rownames(df_norm)
dftemp<-df_norm[,16:19]
dim(dftemp)


dftemp <- dftemp %>%
  filter(gene %in% top_genes$top_genes)

df_temp<-melt(dftemp)

dim(df_temp)
ggplot(df_temp, aes(x = variable, y = gene, color = value)) +
  geom_point() +
  scale_size_continuous(range = c(3, 10)) +  # Tamanho dos pontos
  scale_color_gradient(low = "blue", high = "red") +  # Cor dos pontos de acordo com a expressão
  labs(x = "P-value", y = "Gene", title = "Dotplot de Expressão Gênica vs P-value") +
  theme_minimal() +
  theme(legend.position = "bottom")


dftemp<-dftemp[,-4]


#norm_z<-dftemp
norm_z<-t(apply(dftemp,1,scale))
colnames(norm_z)<-colnames(dftemp)

library(ComplexHeatmap)
library(RColorBrewer)

ht <-Heatmap(
  norm_z,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  column_labels = colnames(norm_z),
  column_names_gp = gpar(fontsize = 16),
  row_names_gp = gpar(fontsize = 20),
  show_heatmap_legend = TRUE,  # Exibe a legenda
  heatmap_legend_param = list(
    title = NULL,  # Sem título para a legenda
    legend_width = unit(6, "cm"),
    legend_height = unit(6, "cm"),  # Ajusta a altura da legenda
    legend_position = "horizontal",  # Coloca a legenda na parte inferior
    legend_direction = "horizontal",
    border = "black"
  ),
  column_title = "REGULATION OF PROGRAMMED CELL DEATH",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  col = colorRampPalette(brewer.pal(9, "Oranges"))(100),  # Coloração "Oranges"
  column_names_rot = 45  # Rotaciona os nomes das colunas em 45 graus (diagonal)
)
draw(ht, heatmap_legend_side = "bottom")


res<-results(dds, contrast=c('condition', 'Tnf','Controle'))
res$negLogp<- (-log10(res$pvalue))
res$genes<- rownames(res)

res<-as.data.frame(res)


res %>%
  filter(genes == 'FGF1')
