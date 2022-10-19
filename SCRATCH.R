quantile(DataTidyACSRISK6::data_tidy_risk6$risk6, c(0.9, 0.95, 0.975, 0.99, 0.995, 1))
quantile(DataTidyACSRISK6::data_tidy_risk6$risk6, rev(1 - c(0.9, 0.95, 0.975, 0.99, 0.995, 1)))

max_val <- quantile(DataTidyACSRISK6::data_tidy_risk6$risk6, 0.975)
min_val <- quantile(DataTidyACSRISK6::data_tidy_risk6$risk6, 0.025)
library(ggplot2)

p <- ggplot(
  DataTidyACSRISK6::data_tidy_risk6,
  aes(y = risk6)
) +
  geom_boxplot()

cowplot::save_plot(
  filename = "p.png",
  p
)

p_wins <- ggplot(
  DataTidyACSRISK6::data_tidy_risk6 |>
    dplyr::mutate(
      risk6 = pmax(pmin(risk6, max_val), min_val)
    ),
  aes(y = risk6)
) +
  geom_boxplot()

cowplot::save_plot(
  filename = "p_wins.png",
  p_wins
)

UtilsDataRSV::view_cols(
  data_raw
)
