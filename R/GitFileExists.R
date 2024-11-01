require("curl")

GitFileExists <- function(repo, path) {
  token <- Sys.getenv("MKNT_READ")
  if (token == "") {
    stop("MKNT_READ environment variable not set")
  }
  url <- sprintf("https://api.github.com/repos/%s/contents/%s", repo, path)
  h <- curl::new_handle(verbose = FALSE)
  curl::handle_setopt(h, customrequest = "GET")
  curl::handle_setheaders(h, Authorization = paste("token", token))
  res <- curl::curl_fetch_memory(url, handle = h)
  
  # Return:
  res[["status_code"]] == "200"
}
