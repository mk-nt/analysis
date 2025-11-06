.ssh <- new.env(parent = emptyenv())

#' Start ssh session with remote server
#' 
#' @returns `SshSession()` is called for its side-effect of beginning a new
#' ssh session and registering this in the package environment for future
#' calls to the remote.
#' @export
SshSession <- function() {
  sess <- .ssh$session
  valid <- FALSE
  
  if (!is.null(sess)) {
    valid <- !inherits(try(ssh_exec_wait(sess, "echo", std_out = NULL),
                           silent = TRUE), "try-error")
  }
  
  if (!valid) {
    if (!is.null(sess)) {
      # Cleanly close the old session to avoid warnings
      try(ssh_disconnect(sess), silent = TRUE)
    }
    sess <- ConnectSSH()
    .ssh$session <- sess
  }
  
  sess
}
