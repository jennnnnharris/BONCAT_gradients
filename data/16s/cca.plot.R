#cca

AF_cca <- cca(AF_1 ~ mix + Condition(block), meta1)

smry <- summary(AF_cca)
smry

#all_cca <- AF_cca
variance <- AF_cca$CA$eig/AF_cca$tot.chi*100

# Extract the model's adjusted R2
RsquareAdj(AF_cca)$adj.r.squared
# Test whether the model is statistically significant

anova.cca(AF_cca, step = 1000, by="term")

# Type 1: Distances among objects reflect their similarities 	
# Type 2: Angles between variables reflect their correlation

ordiplot(AF_cca, scaling = 1, main = "Allele Frequencies partial cca - Scaling 1")
ordiplot(AF_cca, scaling = 2, main = "Allele Frequencies partial cca - Scaling 2")

#custom plot

# percent varience  
perc <- round(100*(summary(AF_cca)$cont$importance[2, 1:2]), 2)
perc

sc_si <- scores(AF_cca, display="sites", choices=c(1,2), scaling=1)
scores(AF_cca, display="sites", choices=c(1,2), scaling=1)


#all
ordiplot(AF_cca, choices=c(1,2), scaling =1, type="none", main="Partial cca Scaling 1", 
         xlab=paste("cca1 (",round(perc[1],1),"% variance explained)"),
         ylab=paste("PC 1 (",round(perc[2],2),"% variance explained)"))
points(sc_si,
       pch=25,
       cex=1,
       #col= c("navy", "orange", "lightblue", "gray")[as.factor(meta1$trt)],
       #bg= c("navy", "orange", "lightblue", "gray")[as.factor(meta1$trt)])
       col= c("black", "gray50", "grey")[as.factor(meta1$year)],
       bg=  c("black", "gray50", "white")[as.factor(meta1$year)])
ordiellipse(sc_si, meta1$trt,  
            kind = "ehull", conf=0.95, label=T, 
            draw = "polygon",
            border = 0,
            lwd=.1,
            col= c("navy", "orange", "lightblue", "gray"),
            alpha = 50)
#dev.off() 