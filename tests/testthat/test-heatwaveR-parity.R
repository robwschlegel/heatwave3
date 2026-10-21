# heatwaveR is the reference for heatwave3's outputs. Compare events, event
# categories and daily categories with heatwaveR, pixel by pixel, on the
# bundled test data, for heatwaves and cold-spells.
#
# Two deliberate differences (heatwaveR bugs, see CLAUDE.md) are excluded:
# - seasons: heatwave3 uses the peak month (Dec-Feb = Summer); heatwaveR's
#   category() means to but uses calendar quarters, with season = "range";
# - series end: heatwaveR joins up to maxGap trailing non-exceeding days at the
#   end of the series into the last event; heatwave3 does not.

skip_if_not_installed("heatwaveR")
# heatwaveR 1.0.0 takes the climatology window's neighbours across the new
# year from the adjacent years, as heatwave3 does; earlier versions wrap
# within each year
skip_if(utils::packageVersion("heatwaveR") < "1.0.0",
        "needs heatwaveR >= 1.0.0 (climatology window at the new year)")

sst_file <- system.file("extdata/sst_test.nc", package = "heatwave3")
cp <- c("1982-01-01", "2011-12-31")

run_both <- function(cold) {
  stem <- tempfile()
  pct <- if (cold) 10 else 90
  detect3(sst_file, name = stem, climatologyPeriod = cp, pctile = pct,
          coldSpells = cold, category = TRUE, daily = "also", quiet = TRUE)
  g <- hw3_read_sst(sst_file, "")
  last_day <- as.Date(max(g$time_days) - 2440588L)
  ref <- lapply(seq_len(g$nlon * g$nlat) - 1L, function(px) {
    ts <- data.frame(t = as.Date(g$time_days - 2440588L),
                     temp = g$sst[px * g$ntime + seq_len(g$ntime)])
    de <- heatwaveR::detect_event(heatwaveR::ts2clm(ts, climatologyPeriod = cp, pctile = pct),
                                  coldSpells = cold)
    list(lon = g$lon[px %/% g$nlat + 1], lat = g$lat[px %% g$nlat + 1], px = px, de = de,
         last_day = last_day,
         cat = heatwaveR::category(de, S = TRUE, climatology = TRUE))
  })
  list(stem = stem, ref = ref, last_day = last_day)
}

# Events of one pixel from both packages, matched by number, without the
# series-end events affected by heatwaveR's trailing-gap join
paired_events <- function(r, ev3) {
  m <- merge(r$de$event, ev3[ev3$pixel_index == r$px, ], by = "event_no",
             suffixes = c(".R", ".3"))
  m[!(m$date_end.R == r$last_day & m$date_end.R != m$date_end.3), ]
}

for (cold in c(FALSE, TRUE)) {
  label <- if (cold) "cold-spells" else "heatwaves"
  runs <- run_both(cold)

  test_that(paste("climatology matches heatwaveR::ts2clm():", label), {
    cd <- hw3_read_clim_nc(paste0(runs$stem, "_clim.nc"))
    g <- hw3_read_sst(sst_file, "")
    for (r in runs$ref) {
      ts <- data.frame(t = as.Date(g$time_days - 2440588L),
                       temp = g$sst[r$px * g$ntime + seq_len(g$ntime)])
      cl <- heatwaveR::ts2clm(ts, climatologyPeriod = cp, pctile = if (cold) 10 else 90)
      cl <- cl[!duplicated(cl$doy), ]
      cl <- cl[order(cl$doy), ]
      expect_identical(cd$seas[r$px * 366 + cl$doy], cl$seas)
      expect_identical(cd$thresh[r$px * 366 + cl$doy], cl$thresh)
    }
  })

  test_that(paste("events match heatwaveR:", label), {
    ev3 <- hw3_export(paste0(runs$stem, "_events.nc"))
    for (r in runs$ref) {
      m <- paired_events(r, ev3)
      expect_gt(nrow(m), 0)
      expect_identical(nrow(r$de$event) - nrow(m) <= 1L, TRUE)  # at most the last event
      for (v in c("date_start", "date_peak", "date_end", "duration"))
        expect_identical(as.numeric(m[[paste0(v, ".3")]]), as.numeric(m[[paste0(v, ".R")]]),
                         label = paste(label, v))
      # With identical climatologies, every metric is identical to heatwaveR's
      for (v in c("intensity_mean", "intensity_max", "intensity_var", "intensity_cumulative",
                  "intensity_mean_relThresh", "intensity_max_relThresh",
                  "intensity_var_relThresh", "intensity_cumulative_relThresh",
                  "intensity_mean_abs", "intensity_max_abs", "intensity_var_abs",
                  "intensity_cumulative_abs", "rate_onset", "rate_decline"))
        expect_identical(as.numeric(m[[paste0(v, ".3")]]), as.numeric(m[[paste0(v, ".R")]]),
                         label = paste(label, v))
    }
  })

  test_that(paste("event categories match heatwaveR::category():", label), {
    inline <- category3(paste0(runs$stem, "_events.nc"))
    # The same events detected without inline categories: computed from the SST
    stem2 <- tempfile()
    file.copy(paste0(runs$stem, "_clim.nc"), paste0(stem2, "_clim.nc"))
    detect_event3(sst_file, name = stem2, coldSpells = cold, quiet = TRUE)
    from_sst <- category3(paste0(stem2, "_events.nc"), paste0(stem2, "_clim.nc"),
                          sst_file = sst_file)
    ev3 <- hw3_export(paste0(runs$stem, "_events.nc"))
    for (c3 in list(inline, from_sst)) {
      for (r in runs$ref) {
        keep <- paired_events(r, ev3)$event_no
        mine <- c3[abs(c3$lon - round(r$lon, 4)) < 1e-6 & abs(c3$lat - round(r$lat, 4)) < 1e-6, ]
        m <- merge(r$cat$event, mine, by = "event_no", suffixes = c(".R", ".3"))
        m <- m[m$event_no %in% keep, ]
        expect_gt(nrow(m), 0)
        expect_identical(m$category.3, m$category.R)
        for (v in c("p_moderate", "p_strong", "p_severe", "p_extreme"))
          expect_identical(m[[paste0(v, ".3")]], m[[paste0(v, ".R")]], label = paste(label, v))
        expect_equal(m$intensity_max, m$i_max, label = paste(label, "i_max"))
      }
    }
  })

  test_that(paste("daily categories match heatwaveR::category(climatology = TRUE):", label), {
    d3 <- hw3_export(paste0(runs$stem, "_events_daily.nc"))
    ev3 <- hw3_export(paste0(runs$stem, "_events.nc"))
    labels <- c("I Moderate", "II Strong", "III Severe", "IV Extreme")
    for (r in runs$ref) {
      clim <- r$cat$climatology
      clim <- clim[clim$event_no %in% paired_events(r, ev3)$event_no, ]
      mine <- d3[d3$lon == r$lon & d3$lat == r$lat & d3$event, ]
      m <- merge(clim, mine, by = "t", suffixes = c(".R", ".3"))
      expect_gt(nrow(m), 0)
      expect_identical(nrow(m), nrow(clim))
      expect_identical(labels[m$category.3], m$category.R)
    }
  })
}
