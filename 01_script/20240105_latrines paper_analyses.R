# This script replicate analyses and figures of the manuscript "Communal tapir latrines are foraging sites for tropical forest vertebrates"
# Submitted to Global Ecology and Conservation Journal on January 05, 2024.
# For more information, contact the corresponding author.
# Lais Lautenschlager - rodrigues.eco12@gmail.com


# loading the packages ----

if(!require(tidyverse))install.packages("tidyverse", dependencies = TRUE)
if(!require(janitor))install.packages("janitor", dependencies = TRUE)
if(!require(hms))install.packages("hms", dependencies = TRUE)
if(!require(ggridges))install.packages("ggridges", dependencies = TRUE)
if(!require(stats))install.packages("stats", dependencies = TRUE)
if(!require(lme4))install.packages("lme4", dependencies = TRUE)
if(!require(remotes))install.packages("remotes", dependencies = TRUE)
if(!require(stringr))install.packages("stringr", dependencies = TRUE)
if(!require(lmtest))install.packages("lmtest", dependencies = TRUE)
if(!require(report))install.packages("report", dependencies = TRUE)
if(!require(circular))install.packages("circular", dependencies = TRUE)
if(!require(lubridate))install.packages("lubridate", dependencies = TRUE)
if(!require(ggpubr))install.packages("ggpubr", dependencies = TRUE)


rm(list = ls())



# ################### DATA ################## ----

# loading the data ----

data_frame <- readr::read_csv("00_data/20231213_dataframe.csv")

spp_groups <- read_csv("00_data/20230811_all_species_list.csv") %>% 
  rename(scientific_name = Species)

# # ######################## Spp interacting with the latrines: figures and models ###################### ----

int_spp <- data_frame |>  
  dplyr::filter(!scientific_name == "rat") |>  
  dplyr::filter(interaction == "non-foraging" | interaction == "foraging" ) |>  
  dplyr::group_by(scientific_name, interaction, latr_id) |>  
  dplyr::summarise(n = n()) |>  
  dplyr::ungroup() |>  
  dplyr::filter(scientific_name %in% c("Aramides_saracura", "Chamaeza_campanisona", "Dasypus_novemcinctus",
                                       "Geotrygon_montana", "Guerlinguetus_brasiliensis", "Metachirus_nudicaudatus",
                                       "Odontophorus_capueira","Tinamus_solitarius", "Turdus_albicollis")) 
int_spp

ggplot2::ggplot(data=int_spp, ggplot2::aes(x= reorder(scientific_name, -n), y=log(n), fill=interaction  )) + 
  ggplot2::geom_boxplot(outlier.shape = NA) +
  ggplot2::scale_fill_manual(values = c("#be014a", "#81684a")) +
  ggplot2::theme_bw() + labs(x = " ", y = "Number of records") +
  ggplot2::guides(fill=ggplot2::guide_legend(title="Interaction status")) + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 75, vjust = 1, hjust=1),
                 axis.title = ggplot2::element_text(size = 20, face = "bold"),
                 axis.text.x.bottom = ggplot2::element_text(size = 18, face = "italic"),
                 axis.text.y = ggplot2::element_text(size = 18),
                 legend.position = c(.85, 0.85),
                 legend.title = ggplot2::element_text(size = 18),
                 legend.text = ggplot2::element_text(size = 18),
                 title = ggplot2::element_text(size = 18)) +
  ylim(0, 7)

#ggsave("02_figures/20231212_latr_int_per_species_bp.svg", w = 25, h = 20, dpi = 300, units = "cm")


# GLMMs --------------------------------------------------------------------

#using negative binomial

int_spp_bin_mdl <- data_frame |>  
  dplyr::filter(!scientific_name == "rat") |>  
  dplyr::filter(interaction == "non-foraging" | interaction == "foraging" ) |> 
  dplyr::select(c(7,18,19,20)) |> 
  dplyr::group_by(scientific_name, interaction, latr_id) |>  
  dplyr::summarise(n = n()) |>  
  dplyr::ungroup() |> 
  dplyr::filter(scientific_name %in% c("Aramides_saracura", "Chamaeza_campanisona", "Dasypus_novemcinctus",
                                       "Geotrygon_montana", "Guerlinguetus_brasiliensis", "Metachirus_nudicaudatus",
                                       "Odontophorus_capueira","Tinamus_solitarius", "Turdus_albicollis"))


# standard error
#se_by_group <- int_spp_bin_mdl %>%
#  group_by(scientific_name, interaction) %>%
#  summarize(StandardError = sd(n) / sqrt(n()))
#se_by_group
#write_csv(int_spp_bin_mdl, "00_data/20231212_int_spp_bin_mdl_table.csv")



# "Aramides_saracura" 
int_nint_spp_mdl.pss.as <- lme4::glmer(n ~ interaction + (1|latr_id), 
                                       data = int_spp_bin_mdl |> 
                                         dplyr::filter(scientific_name == "Aramides_saracura"),
                                       family = poisson)

summary(int_nint_spp_mdl.pss.as)

int_nint_spp_mdl.pss.as.nl <- lme4::glmer(n ~ 1 + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Aramides_saracura"),
                                          family = poisson)

summary(int_nint_spp_mdl.pss.as.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.as.nl, int_nint_spp_mdl.pss.as)



# "Chamaeza_campanisona" 
int_nint_spp_mdl.pss.cc <- lme4::glmer.nb(n ~ interaction + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Chamaeza_campanisona"))

summary(int_nint_spp_mdl.pss.cc)


int_nint_spp_mdl.pss.cc.nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), 
                                             data = int_spp_bin_mdl |> 
                                               dplyr::filter(scientific_name == "Chamaeza_campanisona"))
summary(int_nint_spp_mdl.pss.cc.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.cc.nl, int_nint_spp_mdl.pss.cc)


# "Guerlinguetus_brasiliensis"
int_nint_spp_mdl.pss.gg <- lme4::glmer.nb(n ~ interaction + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Guerlinguetus_brasiliensis"))

summary(int_nint_spp_mdl.pss.gg)

int_nint_spp_mdl.pss.gg.nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), 
                                             data = int_spp_bin_mdl |> 
                                               dplyr::filter(scientific_name == "Guerlinguetus_brasiliensis"))
summary(int_nint_spp_mdl.pss.gg.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.gg.nl, int_nint_spp_mdl.pss.gg)

report::report(int_nint_spp_mdl.pss.gg)

# "Tinamus_solitarius"
int_nint_spp_mdl.pss.ts <- lme4::glmer.nb(n ~ interaction + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Tinamus_solitarius"))

summary(int_nint_spp_mdl.pss.ts)

int_nint_spp_mdl.pss.ts.nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), 
                                             data = int_spp_bin_mdl |> 
                                               dplyr::filter(scientific_name == "Tinamus_solitarius"))
summary(int_nint_spp_mdl.pss.ts.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.ts.nl, int_nint_spp_mdl.pss.ts)


# "Turdus_albicollis"
int_nint_spp_mdl.pss.ta <- lme4::glmer.nb(n ~ interaction + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Turdus_albicollis"))

summary(int_nint_spp_mdl.pss.ta)

int_nint_spp_mdl.pss.ta.nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), 
                                             data = int_spp_bin_mdl |> 
                                               dplyr::filter(scientific_name == "Turdus_albicollis"))
summary(int_nint_spp_mdl.pss.ta.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.ta.nl, int_nint_spp_mdl.pss.ta)



# "Odontophorus_capueira"
int_nint_spp_mdl.pss.oc <- lme4::glmer.nb(n ~ interaction + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Odontophorus_capueira"))

summary(int_nint_spp_mdl.pss.oc)

int_nint_spp_mdl.pss.oc.nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), 
                                             data = int_spp_bin_mdl |> 
                                               dplyr::filter(scientific_name == "Odontophorus_capueira"))
summary(int_nint_spp_mdl.pss.oc.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.oc.nl, int_nint_spp_mdl.pss.oc)



# "Metachirus_nudicaudatus"
int_nint_spp_mdl.pss.mn <- lme4::glmer.nb(n ~ interaction + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Metachirus_nudicaudatus"))

summary(int_nint_spp_mdl.pss.mn)

int_nint_spp_mdl.pss.mn.nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), 
                                             data = int_spp_bin_mdl |> 
                                               dplyr::filter(scientific_name == "Metachirus_nudicaudatus"))
summary(int_nint_spp_mdl.pss.mn.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.mn.nl, int_nint_spp_mdl.pss.mn)



# "Geotrygon_montana"
int_nint_spp_mdl.pss.gm <- lme4::glmer.nb(n ~ interaction + (1|latr_id), 
                                          data = int_spp_bin_mdl |> 
                                            dplyr::filter(scientific_name == "Geotrygon_montana"))

summary(int_nint_spp_mdl.pss.gm)

int_nint_spp_mdl.pss.gm.nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), 
                                             data = int_spp_bin_mdl |> 
                                               dplyr::filter(scientific_name == "Geotrygon_montana"))
summary(int_nint_spp_mdl.pss.gm.nl)

lmtest::lrtest(int_nint_spp_mdl.pss.gm.nl, int_nint_spp_mdl.pss.gm)


# Daytime foraging --------------------------------------------------------


tables_sorted <- data_frame %>% 
  janitor::clean_names() %>% 
  dplyr::select(!c(1,4,5,6,10,15:17,24)) %>% 
  dplyr::mutate(time_video = hms::as_hms(time_video)) %>% 
  dplyr::mutate(scientific_name = stringr::str_replace(scientific_name, " ", ""),
                total_time_sec = as.numeric(total_time_sec),
                time_specific_sec = as.numeric(time_specific_sec),
                time_rate = time_specific_sec/total_time_sec, 
                freq_record = total_time_sec/10) %>%
  dplyr::left_join(spp_groups, by = "scientific_name") %>% 
  dplyr::select(c(1:10,18,19, everything())) %>% 
  janitor::clean_names() |> 
  dplyr::filter(interaction %in% c("foraging", "non-foraging")) 


tables_sorted_time_dt_fg <- tables_sorted %>%
  dplyr::mutate(date_video_mmddyyyy = stringr::str_replace_all(date_video_mmddyyyy, "-", "/"),
                date_video_mmddyyyy = as.POSIXct(paste(tables_sorted$date_video_mmddyyyy, tables_sorted$time_video), 
                                                 format="%Y-%m-%d %H:"),
                date_time = date_video_mmddyyyy) %>%
  dplyr::filter(scientific_name == "Guerlinguetus_brasiliensis" |
                  scientific_name == "Turdus_albicollis" |
                  scientific_name == "Odontophorus_capueira" |
                  scientific_name == "Tinamus_solitarius" |
                  scientific_name == "Chamaeza_campanisona") %>% 
  #dplyr::filter(interaction_latrine == 1) %>% 
  dplyr::mutate(scientific_name = stringr::str_replace(scientific_name, "_", " ")) |> 
  dplyr::distinct() |> 
  tidyr::separate(date_time, c("date", "time"), sep = " ") |> 
  dplyr::mutate(scientific_name = forcats::fct_relevel(scientific_name, c("Guerlinguetus brasiliensis", 
                                                                          "Tinamus solitarius",
                                                                          "Turdus albicollis",
                                                                          "Chamaeza campanisona",
                                                                          "Odontophorus capueira"))) |> 
  dplyr::mutate(time_video = as.POSIXct(time_video),
                interaction = forcats::fct_relevel(interaction, c("non-foraging", "foraging")),
                time_specific_sec = tidyr::replace_na(time_specific_sec, 10))

tables_sorted_time_dt_fg

# Foraging time graph
ggplot2::ggplot(tables_sorted_time_dt_fg, ggplot2::aes(x = time_video, y = scientific_name, fill = interaction, point_color = interaction)) +
  ggridges::geom_density_ridges(color = "black", alpha = 0.7, 
                                jittered_points = TRUE, 
                                quantile_lines=TRUE, 
                                quantiles = 2, 
                                scale = .95, 
                                rel_min_height = .01,
                                point_shape = "|",
                                point_size = 3, 
                                size = 0.25,
                                position = ggridges::position_points_jitter(height = 0)) +
  ggplot2::scale_fill_manual(values = c("#81684a", "#be014a"), 
                             labels = c("Non-foraging", "Foraging"), guide = ggplot2::guide_legend(title = NULL)) +
  ggplot2::scale_discrete_manual("point_color", values = c("#81684a", "#be014a"), guide = "none") +
  xlab("Day time (Hours)") + ylab("") +
  ggridges::theme_ridges(center = TRUE) +
  ggplot2::theme(title = ggplot2::element_text(size = 22),
                 axis.title.x = ggplot2::element_text(face = "bold", size = 20),
                 axis.title.y = ggplot2::element_text(face = "bold", size = 20),
                 axis.text.y = ggplot2::element_text(face="italic", size = 18, color = "black"),
                 axis.text.x = ggplot2::element_text(size=18, angle=0, h = 0.5),
                 legend.text = ggplot2::element_text(face="italic", size = 18),
                 legend.title = ggplot2::element_text(size = 18),
                 legend.direction = "horizontal",
                 legend.position = c(0, 0.98)) +
  ggplot2::scale_x_datetime(date_breaks = "2 hours", date_labels = "%H:%M")

#ggsave("02_figures/20231212_foraging_time_day_graph.svg", w = 35, h = 20, dpi = 600, units = "cm")


# Circular analysis - foraging/non-foraging --------------------------------------------------------

tables_sorted_time_dt_fg_mdl <- tables_sorted_time_dt_fg |> 
  dplyr::select(c(3,10,13,14,21)) |> 
  dplyr::mutate(time = as.numeric(sub(":.*", "", time)))
tables_sorted_time_dt_fg_mdl


# Guerlinguetus brasiliensis
tables_sorted_time_dt_fg_mdl_gb <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Guerlinguetus brasiliensis") |> 
  dplyr::mutate(time_dg = (time/24)*360,
                Hour.Dec = time + (time/60))
tables_sorted_time_dt_fg_mdl_gb


foraging_gb <- as.factor(tables_sorted_time_dt_fg_mdl_gb$interaction)
foraging_gb

angles_gb <- circular::circular(tables_sorted_time_dt_fg_mdl_gb$time_dg)
angles_gb


circ_aov_model_gb <- circular::aov.circular(angles_gb, group = foraging_gb, method = "LRT")
summary(circ_aov_model_gb)
circ_aov_model_gb$mu
circ_aov_model_gb$kappa
circ_aov_model_gb$rho
circ_aov_model_gb$p.value
circ_aov_model_gb$statistic
circ_aov_model_gb$df
circ_aov_model_gb


tables_sorted_time_dt_fg_mdl_gb_for <- tables_sorted_time_dt_fg_mdl_gb |> 
  dplyr::filter(interaction == "foraging")


circ_data_gb_for <-  circular(tables_sorted_time_dt_fg_mdl_gb_for$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_gb_for, bins=24, col = "gray50", prop=1.8, cex = 0.7)


tables_sorted_time_dt_fg_mdl_gb_nfor <- tables_sorted_time_dt_fg_mdl_gb |> 
  dplyr::filter(interaction == "non-foraging")


circ_data_gb_nfor <-  circular(tables_sorted_time_dt_fg_mdl_gb_nfor$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_gb_nfor, bins=24, col = "gray50", prop=1.8, cex = 0.7)



circular_gb <- ggplot2::ggplot(tables_sorted_time_dt_fg_mdl_gb, ggplot2::aes(x = time, fill = interaction)) +
  ggplot2::geom_histogram(breaks = seq(0, 24), width = 2, colour = "grey", alpha=0.8) + 
  ggplot2::coord_polar(start = 0) + ggplot2::theme_minimal() + 
  ggplot2::scale_fill_manual(values = c("#81684a", "#be014a"), 
                             labels = c("Non-foraging", "Foraging"), guide = ggplot2::guide_legend(title = NULL)) +
  ylab("") + ggplot2::ggtitle("Events by Time of day - Guerlinguetus") + 
  ggplot2::scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24)) +
  ggplot2::theme(title = ggplot2::element_text(size = 18),
        axis.title.x = ggplot2::element_text(face = "bold", size = 20),
        axis.title.y = ggplot2::element_text(face = "bold", size = 20),
        axis.text.y = ggplot2::element_text(size = 18, color = "black"),
        axis.text.x = ggplot2::element_text(face="italic", size=18, angle=0, h = 0.5),
        legend.text = ggplot2::element_text(face="italic", size = 18),
        legend.title = ggplot2::element_text(size = 16),
        legend.position = c(0.55,0.7))
circular_gb


# Tinamus solitarius

tables_sorted_time_dt_fg_mdl_ts <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Tinamus solitarius") |> 
  dplyr::mutate(time_dg = (time/24)*360,
                Hour.Dec = time + (time/60))
tables_sorted_time_dt_fg_mdl_ts

hist(tables_sorted_time_dt_fg_mdl_ts$time_dg)

foraging_ts <- as.factor(tables_sorted_time_dt_fg_mdl_ts$interaction)
foraging_ts
angles_ts <- circular::circular(tables_sorted_time_dt_fg_mdl_ts$time_dg)
angles_ts


circ_aov_model_ts <- circular::aov.circular(angles_ts, group = foraging_ts, method = "LRT")
summary(circ_aov_model_ts)
circ_aov_model_ts$mu
circ_aov_model_ts$kappa
circ_aov_model_ts$rho
circ_aov_model_ts$p.value
circ_aov_model_ts$statistic
circ_aov_model_ts$df
circ_aov_model_ts


tables_sorted_time_dt_fg_mdl_ts_for <- tables_sorted_time_dt_fg_mdl_ts |> 
  dplyr::filter(interaction == "foraging")


circ_data_ts_for <-  circular(tables_sorted_time_dt_fg_mdl_ts_for$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_ts_for, bins=24, col = "gray50", prop=1.8, cex = 0.7)


tables_sorted_time_dt_fg_mdl_ts_nfor <- tables_sorted_time_dt_fg_mdl_ts |> 
  dplyr::filter(interaction == "non-foraging")


circ_data_ts_nfor <-  circular(tables_sorted_time_dt_fg_mdl_ts_nfor$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_ts_nfor, bins=24, col = "gray50", prop=1.8, cex = 0.7)



circular_ts <- ggplot2::ggplot(tables_sorted_time_dt_fg_mdl_ts, ggplot2::aes(x = time, fill = interaction)) +
  ggplot2::geom_histogram(breaks = seq(0, 24), width = 2, colour = "grey", alpha=0.8) + 
  ggplot2::coord_polar(start = 0) + ggplot2::theme_minimal() + 
  ggplot2::scale_fill_manual(values = c("#81684a", "#be014a"), 
                             labels = c("Non-foraging", "Foraging"), guide = ggplot2::guide_legend(title = NULL)) +
  ylab("") + ggplot2::ggtitle("Events by Time of day - Tinamus") + 
  ggplot2::scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24)) +
  ggplot2::theme(title = ggplot2::element_text(size = 18),
                 axis.title.x = ggplot2::element_text(face = "bold", size = 20),
                 axis.title.y = ggplot2::element_text(face = "bold", size = 20),
                 axis.text.y = ggplot2::element_text(size = 18, color = "black"),
                 axis.text.x = ggplot2::element_text(face="italic", size=18, angle=0, h = 0.5),
                 legend.text = ggplot2::element_text(face="italic", size = 18),
                 legend.title = ggplot2::element_text(size = 16),
                 legend.position = c(0.55,0.7))
circular_ts


# Turdus albicollis


tables_sorted_time_dt_fg_mdl_ta <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Turdus albicollis") |> 
  dplyr::mutate(time_dg = (time/24)*360,
                Hour.Dec = time + (time/60))
tables_sorted_time_dt_fg_mdl_ta


foraging_ta <- as.factor(tables_sorted_time_dt_fg_mdl_ta$interaction)
foraging_ta
angles_ta <- circular::circular(tables_sorted_time_dt_fg_mdl_ta$time_dg)
angles_ta


circ_aov_model_ta <- circular::aov.circular(angles_ta, group = foraging_ta, method = "LRT")
summary(circ_aov_model_ta)
circ_aov_model_ta$mu
circ_aov_model_ta$kappa
circ_aov_model_ta$rho
circ_aov_model_ta$p.value
circ_aov_model_ta$statistic
circ_aov_model_ta$df
circ_aov_model_ta


tables_sorted_time_dt_fg_mdl_ta_for <- tables_sorted_time_dt_fg_mdl_ta |> 
  dplyr::filter(interaction == "foraging")


circ_data_ta_for <-  circular(tables_sorted_time_dt_fg_mdl_ta_for$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_ta_for, bins=24, col = "gray50", prop=1.8, cex = 0.7)


tables_sorted_time_dt_fg_mdl_ta_nfor <- tables_sorted_time_dt_fg_mdl_ta |> 
  dplyr::filter(interaction == "non-foraging")


circ_data_ta_nfor <-  circular(tables_sorted_time_dt_fg_mdl_ta_nfor$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_ta_nfor, bins=24, col = "gray50", prop=1.8, cex = 0.7)


circular_ta <- ggplot2::ggplot(tables_sorted_time_dt_fg_mdl_ta, ggplot2::aes(x = time, fill = interaction)) +
  ggplot2::geom_histogram(breaks = seq(0, 24), width = 2, colour = "grey", alpha=0.8) + 
  ggplot2::coord_polar(start = 0) + ggplot2::theme_minimal() + 
  ggplot2::scale_fill_manual(values = c("#81684a", "#be014a"), 
                             labels = c("Non-foraging", "Foraging"), guide = ggplot2::guide_legend(title = NULL)) +
  ylab("") + ggplot2::ggtitle("Events by Time of day - Turdus") + 
  ggplot2::scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24)) +
  ggplot2::theme(title = ggplot2::element_text(size = 18),
                 axis.title.x = ggplot2::element_text(face = "bold", size = 20),
                 axis.title.y = ggplot2::element_text(face = "bold", size = 20),
                 axis.text.y = ggplot2::element_text(size = 18, color = "black"),
                 axis.text.x = ggplot2::element_text(face="italic", size=18, angle=0, h = 0.5),
                 legend.text = ggplot2::element_text(face="italic", size = 18),
                 legend.title = ggplot2::element_text(size = 16),
                 legend.position = c(0.55,0.7))
circular_ta


# Chamaeza campanisona

tables_sorted_time_dt_fg_mdl_cc <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Chamaeza campanisona") |> 
  dplyr::mutate(time_dg = (time/24)*360,
                Hour.Dec = time + (time/60))
tables_sorted_time_dt_fg_mdl_cc


foraging_cc <- as.factor(tables_sorted_time_dt_fg_mdl_cc$interaction)
foraging_cc
angles_cc <- circular::circular(tables_sorted_time_dt_fg_mdl_cc$time_dg)
angles_cc


circ_aov_model_cc <- circular::aov.circular(angles_cc, group = foraging_cc, method = "LRT")
summary(circ_aov_model_cc)
circ_aov_model_cc$mu
circ_aov_model_cc$kappa
circ_aov_model_cc$rho
circ_aov_model_cc$p.value
circ_aov_model_cc$statistic
circ_aov_model_cc$df
circ_aov_model_cc



tables_sorted_time_dt_fg_mdl_cc_for <- tables_sorted_time_dt_fg_mdl_cc |> 
  dplyr::filter(interaction == "foraging")


circ_data_cc_for <-  circular(tables_sorted_time_dt_fg_mdl_cc_for$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_cc_for, bins=24, col = "gray50", prop=1.8, cex = 0.7)


tables_sorted_time_dt_fg_mdl_cc_nfor <- tables_sorted_time_dt_fg_mdl_cc |> 
  dplyr::filter(interaction == "non-foraging")


circ_data_cc_nfor <-  circular(tables_sorted_time_dt_fg_mdl_cc_nfor$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_cc_nfor, bins=24, col = "gray50", prop=1.8, cex = 0.7)


circular_cc <- ggplot2::ggplot(tables_sorted_time_dt_fg_mdl_cc, ggplot2::aes(x = time, fill = interaction)) +
  ggplot2::geom_histogram(breaks = seq(0, 24), width = 2, colour = "grey", alpha=0.8) + 
  ggplot2::coord_polar(start = 0) + ggplot2::theme_minimal() + 
  ggplot2::scale_fill_manual(values = c("#81684a", "#be014a"), 
                             labels = c("Non-foraging", "Foraging"), guide = ggplot2::guide_legend(title = NULL)) +
  ylab("") + ggplot2::ggtitle("Events by Time of day - Chamaeza") + 
  ggplot2::scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24)) +
  ggplot2::theme(title = ggplot2::element_text(size = 18),
                 axis.title.x = ggplot2::element_text(face = "bold", size = 20),
                 axis.title.y = ggplot2::element_text(face = "bold", size = 20),
                 axis.text.y = ggplot2::element_text(size = 18, color = "black"),
                 axis.text.x = ggplot2::element_text(face="italic", size=18, angle=0, h = 0.5),
                 legend.text = ggplot2::element_text(face="italic", size = 18),
                 legend.title = ggplot2::element_text(size = 16),
                 legend.position = c(0.55,0.7))
circular_cc




# Odontophorus capueira

tables_sorted_time_dt_fg_mdl_oc <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Odontophorus capueira") |> 
  dplyr::mutate(time_dg = (time/24)*360,
                Hour.Dec = time + (time/60))
tables_sorted_time_dt_fg_mdl_oc


foraging_oc <- as.factor(tables_sorted_time_dt_fg_mdl_oc$interaction)
foraging_oc
angles_oc <- circular::circular(tables_sorted_time_dt_fg_mdl_oc$time_dg)
angles_oc


circ_aov_model_oc <- circular::aov.circular(angles_oc, group = foraging_oc, method = "LRT")
summary(circ_aov_model_oc)
circ_aov_model_oc$mu
circ_aov_model_oc$kappa
circ_aov_model_oc$rho
circ_aov_model_oc$p.value
circ_aov_model_oc$statistic
circ_aov_model_oc$df
circ_aov_model_oc



tables_sorted_time_dt_fg_mdl_oc_for <- tables_sorted_time_dt_fg_mdl_oc |> 
  dplyr::filter(interaction == "foraging")


circ_data_oc_for <-  circular(tables_sorted_time_dt_fg_mdl_oc_for$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_oc_for, bins=24, col = "gray50", prop=1.8, cex = 0.7)


tables_sorted_time_dt_fg_mdl_oc_nfor <- tables_sorted_time_dt_fg_mdl_oc |> 
  dplyr::filter(interaction == "non-foraging")


circ_data_oc_nfor <-  circular(tables_sorted_time_dt_fg_mdl_oc_nfor$Hour.Dec, type="angles",units="hours",template="clock24")
rose.diag(circ_data_oc_nfor, bins=24, col = "gray50", prop=1.8, cex = 0.7)


circular_oc <- ggplot2::ggplot(tables_sorted_time_dt_fg_mdl_oc, ggplot2::aes(x = Hour.Dec, fill = interaction)) +
  ggplot2::geom_histogram(breaks = seq(0, 24), width = 2, colour = "grey", alpha=0.8) + 
  ggplot2::coord_polar(start = 0) + ggplot2::theme_minimal() + 
  ggplot2::scale_fill_manual(values = c("#81684a", "#be014a"), 
                             labels = c("Non-foraging", "Foraging"), guide = ggplot2::guide_legend(title = NULL)) +
  ylab("") + ggplot2::ggtitle("Events by Time of day - Odontophorus") + 
  ggplot2::scale_x_continuous("", limits = c(0, 24), breaks = seq(0, 24), labels = seq(0, 24)) +
  ggplot2::theme(title = ggplot2::element_text(size = 18),
                 axis.title.x = ggplot2::element_text(face = "bold", size = 20),
                 axis.title.y = ggplot2::element_text(face = "bold", size = 20),
                 axis.text.y = ggplot2::element_text(size = 18, color = "black"),
                 axis.text.x = ggplot2::element_text(face="italic", size=18, angle=0, h = 0.5),
                 legend.text = ggplot2::element_text(face="italic", size = 18),
                 legend.title = ggplot2::element_text(size = 16),
                 legend.position = c(0.55,0.7))
circular_oc

# merging the graphs

figure <- ggpubr::ggarrange(circular_gb, circular_ta, circular_ts, circular_cc, circular_oc,
                    ncol = 5, nrow = 1)
figure

#ggplot2::ggsave("02_figures/circular_plots.svg", w = 60, h = 40, units = "cm", dpi = 500)



########################### organizing the data for foraging after defecation events ##########################

# linear regression graph of interactions after pooping by counting number  -------------------------------------

events_int <- readr::read_csv("00_data/20231213_interactions_and_non_interactions_after_defecation.csv")


# names(events_int)
# 
# poopinteing_events_int <- events_int |> 
#   dplyr::filter(Days_after_defecation <= 30) |> 
#   dplyr::mutate(Species = stringr::str_replace(Species, "_", " ")) |> 
#   dplyr::distinct()
# poopinteing_events_int
# names(poopinteing_events_int)
# 
# 
# ggplot2::ggplot(poopinteing_events_int, ggplot2::aes(x = Days_after_defecation,  y = log(Interactions_Seconds), 
#                                                  color = Interaction_type, 
#                                                  fill = Interaction_type)) +
#   ggplot2::geom_point() +
#   ggplot2::geom_smooth(method = "loess", se=TRUE) +
#   ggplot2::scale_fill_manual(values = c("#be014a","#81684a"), 
#                              labels = c("Foraging","Non-foraging"), guide = ggplot2::guide_legend(title = NULL)) +
#   ggplot2::scale_color_manual(values = c("#be014a","#81684a"), 
#                               labels = c("Foraging","Non-foraging"), guide = ggplot2::guide_legend(title = NULL)) +
#   #geom_rug() +
#   xlab("Days") + ylab("Number of records (log)") +
#   ggplot2::theme_bw() +
#   ggplot2::theme(axis.title.x = ggplot2::element_text(face = "bold", size = 18),
#                  axis.title.y = ggplot2::element_text(face = "bold", size = 18),
#                  axis.text.y = ggplot2::element_text(size = 15),
#                  axis.text.x = ggplot2::element_text(size=15, angle=0, h = 0.5),
#                  legend.text = ggplot2::element_text(face="italic", size = 18),
#                  strip.text = ggplot2::element_text(size=12, face="italic"),
#                  legend.position = c(0.7, 0.2),
#   ) +
#   ggplot2::scale_x_continuous(breaks=seq(0, 30, 5)) +
#   ggplot2::facet_wrap(.~Species, nrow = 3, ncol = 2, scales = "free")
# 
# #ggsave("02_figures/20231118_pooping_events_int_opt6.svg", w = 15, h = 20, dpi = 400, units = "cm")

# linear regression graph of interactions after pooping by binomical data  -------------------------------------

poopinteing_events_int_bin <- events_int |> 
  dplyr::filter(Days_after_defecation <= 30) |> 
  dplyr::mutate(Species = stringr::str_replace(Species, "_", " "),
                Interactions_foraging = if_else(Interaction_type == "foraging", 1, 0),
                Interactions_nonforaging = if_else(Interaction_type == "non-foraging", 1, 0),
                Date_Interaction = as.factor(Date_Interaction)) 

names(poopinteing_events_int_bin)


# "Guerlinguetus brasiliensis"

poopinteing_events_int_bin_gg <- poopinteing_events_int_bin |> 
  dplyr::filter(Species == "Guerlinguetus brasiliensis")

bn_model_gg <- lme4::glmer(Interactions_foraging ~ Days_after_defecation +  (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_gg )
summary(bn_model_gg)

bn_model_gg_nl <- lme4::glmer(Interactions_foraging ~ 1 +  (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_gg)
summary(bn_model_gg_nl)

lmtest::lrtest(bn_model_gg_nl, bn_model_gg)


poopinteing_events_int_bin_gg$predicted_probs <- predict(bn_model_gg, type = "response")
poopinteing_events_int_bin_gg


#"Tinamus solitarius"

poopinteing_events_int_bin_ts <- poopinteing_events_int_bin |> 
  dplyr::filter(Species == "Tinamus solitarius")

bn_model_ts <- lme4::glmer(Interactions_foraging ~ Days_after_defecation + (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_ts)
summary(bn_model_ts)

bn_model_ta_nl <- lme4::glmer(Interactions_foraging ~ 1 +  (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_ts)
summary(bn_model_ta_nl)

lmtest::lrtest(bn_model_ta_nl, bn_model_ts)

poopinteing_events_int_bin_ts$predicted_probs <- predict(bn_model_ts, type = "response")
poopinteing_events_int_bin_ts


# "Turdus albicollis"

poopinteing_events_int_bin_ta <- poopinteing_events_int_bin |> 
  dplyr::filter(Species == "Turdus albicollis")

bn_model_ta <- lme4::glmer(Interactions_foraging ~ Days_after_defecation +  (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_ta)
summary(bn_model_ta)
#report(bn_model_ta)

bn_model_ta_nl <- lme4::glmer(Interactions_foraging ~ 1 +  (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_ta)
summary(bn_model_ta_nl)

lmtest::lrtest(bn_model_ta_nl, bn_model_ta)

poopinteing_events_int_bin_ta$predicted_probs <- predict(bn_model_ta, type = "response")
poopinteing_events_int_bin_ta


# "Chamaeza campanisona"

poopinteing_events_int_bin_cc <- poopinteing_events_int_bin |> 
  dplyr::filter(Species == "Chamaeza campanisona")

bn_model_cc <- lme4::glmer(Interactions_foraging ~ Days_after_defecation + (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_cc)
summary(bn_model_cc)
#report(bn_model_cc)

bn_model_cc_nl <- lme4::glmer(Interactions_foraging ~ 1 + (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_cc)
summary(bn_model_cc_nl)

lmtest::lrtest(bn_model_cc_nl, bn_model_cc)

poopinteing_events_int_bin_cc$predicted_probs <- predict(bn_model_cc, type = "response")
poopinteing_events_int_bin_cc


# "Odontophorus capueira"

poopinteing_events_int_bin_oc <- poopinteing_events_int_bin |> 
  dplyr::filter(Species == "Odontophorus capueira")

bn_model_oc <- lme4::glmer(Interactions_foraging ~ Days_after_defecation + (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_oc)
summary(bn_model_oc)

bn_model_oc_nl <- lme4::glmer(Interactions_foraging ~ 1 +  (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_oc)
summary(bn_model_oc_nl)

lmtest::lrtest(bn_model_oc_nl, bn_model_oc)

poopinteing_events_int_bin_oc$predicted_probs <- predict(bn_model_oc, type = "response")
poopinteing_events_int_bin_oc


# Graph binomial regression

poopinteing_events_int_bin_gg_predictors <- rbind(poopinteing_events_int_bin_gg, poopinteing_events_int_bin_ts, poopinteing_events_int_bin_ta, poopinteing_events_int_bin_cc, poopinteing_events_int_bin_oc)
names(poopinteing_events_int_bin_gg_predictors)

ggplot(poopinteing_events_int_bin_gg_predictors, aes(x = Days_after_defecation, y = predicted_probs)) + 
  geom_point(alpha=.5, size = 2, color="black") + 
  stat_smooth(method="glm", se=TRUE, method.args = list(family=binomial), color = "#4682B4") +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=14),
        strip.text.x = element_text(size = 12, face = "italic")) +
  ylab("Foraging probability") + xlab("Days") +
  #scale_y_continuous(breaks = seq(0, 1)) +
  facet_wrap(.~Species, ncol = 5) +
  ggplot2::scale_x_continuous(breaks=seq(0, 30, 5))


#ggsave("02_figures/20231213_regression_line.png", w = 35, h = 10, dpi = 500, units = "cm")

# end ------------------------------------------------------------------------------------