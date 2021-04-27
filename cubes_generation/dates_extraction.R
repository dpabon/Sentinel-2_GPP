dates_extraction <- function(path, pattern){
  require(stringr)
  require(lubridate)
  dates <- list.files(path = path, pattern = glob2rx(pattern))
  dates <- str_sub(dates, start = 1L, end = 8L)
  dates <- as.Date(dates, "%Y%m%d")
  return(dates)
}
