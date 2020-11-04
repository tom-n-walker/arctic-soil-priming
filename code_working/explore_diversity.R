rm(list = ls())
library(drake)

loadd()

met_nums = metabolites %>% 
  select(-sampleID)

richness = data.frame(sampleID = metabolites$sampleID,
                      richness = rowSums(met_nums > 0),
                      shannon = diversity(met_nums),
                      simpson = diversity(met_nums, "simpson"))

all = left_join(priming, richness)
all = all[all$Rprotein < 4, ]

m = lm(richness ~ horizon, all)
anova(m)
contrast(emmeans(m, ~ horizon), "pairwise")

sig = data.frame(horizon = c("O", "A", "B", "F", "J"),
                 signif = c("a", "a", "b", "b", "a"))
sig$horizon = factor(sig$horizon, levels = c("O", "A", "B", "F", "J"))
all$horizon = factor(all$horizon, levels = c("O", "A", "B", "F", "J"))


plot_horizons = ggplot(all) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = horizon, y = richness, fill = horizon) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_text(data = sig, aes(label = signif, y = 150)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.shape = NA) +
  guides(fill = "none") +
  xlab("Soil horizon") + ylab("Soil organic matter richness")

plot_priming = ggplot(all) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = richness, y = Rprotein * 100, fill = horizon, col = horizon) +
  geom_smooth(method = "lm", se = F, size = 0.25) +
  geom_smooth(aes(col = NULL, fill = NULL), method = "lm", se = F, size = 0.5, col = "black") +
  geom_point(shape = 21, col = "black") +
  scale_color_brewer(palette = "Pastel1") +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none", col = "none") +
  xlab("Soil organic matter richness") +
  ylab("Priming effect (%)")

postscript(file = "om_diversity_priming.eps", width = mm2in(180), height = mm2in(90))
cowplot::plot_grid(plot_horizons, plot_priming)
dev.off()

plot(Rprotein ~ richness, all)

anova(lm(Rprotein ~ horizon * richness, all))


no_na = all[!is.na(all$richness), ]


anova(lm(Rprotein ~ poly(richness, 2), no_na))

ggplot(all) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  aes(x = richness, y = Rprotein) +
  geom_point() +
