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
# 
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


# correlations

tables_sorted_time_dt_fg_mdl <- tables_sorted_time_dt_fg |> 
  dplyr::select(c(3,10,13,14,21)) |> 
  dplyr::mutate(time = as.numeric(sub(":.*", "", time))) |> 
  ggpubr::group_by(latr_id, scientific_name, interaction, time) |> 
  dplyr::summarise(n = n()) |> 
  dplyr::ungroup()
tables_sorted_time_dt_fg_mdl


# Guerlinguetus brasiliensis
tables_sorted_time_dt_fg_mdl_gb <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Guerlinguetus brasiliensis")
tables_sorted_time_dt_fg_mdl_gb


tables_sorted_time_dt_fg_mdl_gb_ng <- lme4::glmer.nb(n ~ interaction * time + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_gb)
summary(tables_sorted_time_dt_fg_mdl_gb_ng)

tables_sorted_time_dt_fg_mdl_gb_ng_nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_gb)
summary(tables_sorted_time_dt_fg_mdl_gb_ng_nl)

lmtest::lrtest(tables_sorted_time_dt_fg_mdl_gb_ng_nl, tables_sorted_time_dt_fg_mdl_gb_ng)



# Tinamus solitarius
tables_sorted_time_dt_fg_mdl_ts <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Tinamus solitarius")
tables_sorted_time_dt_fg_mdl_ts

tables_sorted_time_dt_fg_mdl_ts_ng <- lme4::glmer.nb(n ~ interaction * time + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_ts)
summary(tables_sorted_time_dt_fg_mdl_ts_ng)

tables_sorted_time_dt_fg_mdl_ts_ng_nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_ts)
summary(tables_sorted_time_dt_fg_mdl_ts_ng_nl)

lmtest::lrtest(tables_sorted_time_dt_fg_mdl_ts_ng_nl, tables_sorted_time_dt_fg_mdl_ts_ng)



# Turdus albicollis
tables_sorted_time_dt_fg_mdl_ta <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Turdus albicollis")
tables_sorted_time_dt_fg_mdl_ta


tables_sorted_time_dt_fg_mdl_ta_ng <- lme4::glmer.nb(n ~ interaction * time + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_ta)
summary(tables_sorted_time_dt_fg_mdl_ta_ng)

tables_sorted_time_dt_fg_mdl_ta_ng_nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_ta)
summary(tables_sorted_time_dt_fg_mdl_ta_ng_nl)

lmtest::lrtest(tables_sorted_time_dt_fg_mdl_ta_ng_nl, tables_sorted_time_dt_fg_mdl_ta_ng)



# Chamaeza campanisona

tables_sorted_time_dt_fg_mdl_cc <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Chamaeza campanisona")
tables_sorted_time_dt_fg_mdl_cc

tables_sorted_time_dt_fg_mdl_cc_ng <- lme4::glmer.nb(n ~ interaction * time + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_cc)
summary(tables_sorted_time_dt_fg_mdl_cc_ng)

tables_sorted_time_dt_fg_mdl_cc_ng_nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_cc)
summary(tables_sorted_time_dt_fg_mdl_cc_ng_nl)

lmtest::lrtest(tables_sorted_time_dt_fg_mdl_cc_ng_nl, tables_sorted_time_dt_fg_mdl_cc_ng)



# Odontophorus capueira

tables_sorted_time_dt_fg_mdl_oc <- tables_sorted_time_dt_fg_mdl |> 
  dplyr::filter(scientific_name == "Odontophorus capueira")
tables_sorted_time_dt_fg_mdl_oc


tables_sorted_time_dt_fg_mdl_oc_ng <- lme4::glmer.nb(n ~ interaction * time + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_oc)
summary(tables_sorted_time_dt_fg_mdl_oc_ng)

tables_sorted_time_dt_fg_mdl_oc_ng_nl <- lme4::glmer.nb(n ~ 1 + (1|latr_id), data = tables_sorted_time_dt_fg_mdl_oc)
summary(tables_sorted_time_dt_fg_mdl_oc_ng_nl)

lmtest::lrtest(tables_sorted_time_dt_fg_mdl_oc_ng_nl, tables_sorted_time_dt_fg_mdl_oc_ng)



########################### organizing the data for poop decay ##########################




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
# 
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

bn_model_ts_nl <- lme4::glmer(Interactions_foraging ~ 1 +  (1|LATR_ID/Date_Interaction), family = "binomial", data = poopinteing_events_int_bin_ts)
summary(bn_model_ts_nl)

lmtest::lrtest(bn_model_ts_nl, bn_model_ts)

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

