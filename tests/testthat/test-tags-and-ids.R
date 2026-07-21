test_that("importTags harmonises tag metadata and parses tagging dates", {
  tg <- importTags(data.frame(
    Transmitter = c("A69-1602-111", "A69-1602-222"), ID = c("shark_01", "shark_02"),
    Tagdeployed = c("2023-06-01", "2023-06-02"), Sex = c("F", "M"), Tl_cm = c(120, 135)),
    source = "vue")
  expect_true(all(c("ID", "transmitter", "tagging_date") %in% colnames(tg)))
  expect_s3_class(tg$tagging_date, "POSIXct")
  expect_true(is.numeric(tg$length))
})

test_that("importTags combines GLATOS codespace + id", {
  tg <- importTags(data.frame(animal_id = "f1", transmitter_codespace = "A69-1602",
                              transmitter_id = "111", utc_release_date_time = "2023-06-01",
                              tag_serial_number = "1450101"), source = "glatos")
  expect_equal(tg$transmitter, "A69-1602-111")
  expect_equal(tg$ID, "f1")
})

test_that("assignAnimalIDs maps transmitters to animal IDs and attaches tagging dates", {
  det <- as_moby(data.frame(
    ID = factor(c("A69-1602-111", "A69-1602-222")),
    transmitter = c("A69-1602-111", "A69-1602-222"),
    datetime = as.POSIXct(c("2023-06-05", "2023-06-06"), tz = "UTC"),
    station = "R1", lon = -8, lat = 37), epsg.code = 32629)
  tags <- importTags(data.frame(Transmitter = c("A69-1602-111", "A69-1602-222"),
                                ID = c("shark_01", "shark_02"),
                                Tagdeployed = c("2023-06-01", "2023-06-02"),
                                Sex = c("F", "M"), Tl_cm = c(120, 135)), source = "vue")
  res <- assignAnimalIDs(det, tags, keep.cols = c("sex", "length"))
  expect_true(all(as.character(res$ID) %in% c("shark_01", "shark_02")))
  expect_true(all(c("sex", "length") %in% colnames(res)))
  expect_false(is.null(mobyMeta(res)$tagging.dates))
  expect_equal(mobyMeta(res)$epsg.code, 32629)  # original metadata preserved
})

test_that("assignAnimalIDs matches by trailing numeric code and warns on unmatched", {
  det <- as_moby(data.frame(ID = factor(c("A69-1602-111", "A69-1602-999")),
                            transmitter = c("A69-1602-111", "A69-1602-999"),
                            datetime = as.POSIXct(c("2023-06-05", "2023-06-06"), tz = "UTC")))
  tags <- importTags(data.frame(Transmitter = "111", ID = "s1", Tagdeployed = "2023-06-01"),
                     source = "vue")
  res <- suppressWarnings(assignAnimalIDs(det, tags))
  expect_equal(as.character(res$ID[res$transmitter == "A69-1602-111"]), "s1")  # numeric match
  expect_true(any(is.na(res$ID)))                                              # 999 unmatched
  expect_warning(assignAnimalIDs(det, tags), "not found")
})
