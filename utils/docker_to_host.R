# docker_to_host.R
#
# Purpose: Helper functions to communicate with Cytoscape EM from R in docker
# Version: 1.0
# Date: 2021-04-01
# Author: Feifei Li <ff.li@mail.utoronto.ca>
#
# Dependencies: BiocManager 3.12 or higher, httr, RCy3, htmltools
#
# ToDo: Test docker_to_unix, docker_to_em, get_em_unix under UNIX environment.
# Notes:
#     docker_to_dos is for Windows only. docker_to_unix is for UNIX OS only.
#     docker_to_em requires EM 3.3.2 or higher version.
#
# ==============================================================================


# ====  PACKAGES  ==============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
    BiocManager::install(version = "3.12", ask = FALSE)
}

# Required for docker_to_em to communicate with CyREST
if (! requireNamespace("httr", quietly=TRUE))
    BiocManager::install("httr", ask = FALSE)

# Required to get screenshot of networks from Cytoscape
if (! requireNamespace("RCy3", quietly = TRUE))
    BiocManager::install("RCy3", ask = FALSE)

# Required to embed image files in html notebook
if (! requireNamespace("htmltools", quietly = TRUE))
    install.packages("htmltools")

# ====  FUNCTIONS  =============================================================

docker_to_dos <- function(localPath, host_dir) {
	# Purpose:
	#     Given an absolute path to a file in docker,
    #     map its path to the local host machine w/ DOS file system.
    #     Intended for Windows users
	# Parameters:
	#     localPath: absolute path to a file in docker
	#     host_dir: absolute path to the directory on the host loacl machine
	# Value:
	#     result: hostPath
    
    docker_dir   <- getwd()
	
    hostPath     <- gsub(getwd(), host_dir, localPath)
    hostPath     <- gsub(pattern     = "/+",
                         replacement = "\\\\",
                         hostPath)

	return(hostPath)
}


docker_to_unix <- function(localPath, host_dir) {
  # Purpose:
  #     Given an absolute path to a local file in docker,
  #     directly upload it to EM app in Cytoscape.
  #     For UNIX/UNIX-like OS only.
  # Parameters:
  #     localPath: absolute path to a file in docker
  #     host_path: absolute path to the directory on the host loacl machine
  # Value:
  #     result: hostPath
  
  docker_dir <- dir(getwd())
  
  host_dir   <- gsub(pattern     = docker_dir,
                     replacement = host_dir,
                     localPath)
  
  return(hostPath)
}


docker_to_em <- function(localPath, port = 1234) {
    # Purpose:
    #     Communicate with Cytoscape directly via CyREST.
    #     No need to provide local host path.
    #     For UNIX/UNIX-like OS only. Required EM 3.3.2 Beta version.
    # Parameters:
    #     localPath: absolute path to a file in docker
    #     port: host port to access EM in cytoscape, default 1234.
    # Value:
    #     result: hostPath
    
    URL      <- paste0("http://localhost:", port,
                       "/enrichmentmap/textfileupload?fileName=")
    fname    <- basename(localPath)
    request  <- httr::POST(url    = paste0(URL, fname),
                           config = list(),
                           body   = list(file = httr::upload_file(localPath)),
                           encode = "multipart")
    hostPath <- content(request, "parsed")$path
    
    return(hostPath)
}


get_em_win <- function(title, host_dir, network = NULL) {
    # Purpose:
    #     Export a screenshot of the current network from docker to Windows host
    # Parameters:
    #     title: Title of the network
    #     host_dir: working directory of the host machine
    #     network: Cytoscape suid of a network
    # Value:
    #     network_screenshot: screen shot of the EM
    
    if (!dir.exists(file.path(getwd(), "screenshots")))
        dir.create(file.path(getwd(), "screenshots"))
    
    fname <- gsub("[^a-zA-Z0-9]", "_", title)
    fname <- paste0(fname, ".png")
    
    export_network_img_docker <- file.path(getwd(), "screenshots", fname)
    export_network_img_host   <- docker_to_dos(export_network_img_docker,
                                               host_dir)
    
    if (file.exists(export_network_img_docker))
        response <- file.remove(export_network_img_docker)
    
    suid     <- if (is.null(network)) NULL else as.numeric(network)
    response <- RCy3::exportImage(export_network_img_host,
                                  network = suid,
                                  type = "png")
    
    network_screenshot <- htmltools::img(
        src   = knitr::image_uri(export_network_img_docker),
        alt   = title,
        style = "margin:0px auto;display:block"
    )
    
    return(network_screenshot)
}


get_em_unix <- function(title, host_dir, network = NULL) {
    # Purpose:
    #     Export a screenshot of the current network from docker to UNIX host
    # Parameters:
    #     title: Title of the network
    #     host_dir: working directory of the host machine
    #     network: Cytoscape suid of a network
    # Value:
    #     network_screenshot: screen shot of the EM
    
    if (!dir.exists(file.path(getwd(), "screenshots")))
        dir.create(file.path(getwd(), "screenshots"))
    
    fname <- gsub("[^a-zA-Z0-9]", "_", title)
    fname <- paste0(fname, ".png")
    
    export_network_img_docker <- file.path(getwd(), "screenshots", fname)
    export_network_img_host   <- gsub(pattern     = getwd(),
                                      replacement = host_dir,
                                      export_network_img_docker)
    
    if (file.exists(export_network_img_docker))
        response <- file.remove(export_network_img_docker)
    
    suid     <- if (is.null(network)) NULL else as.numeric(network)
    response <- RCy3::exportImage(export_network_img_host,
                                  network = suid,
                                  type = "png")
    
    network_screenshot <- htmltools::img(
        src   = knitr::image_uri(export_network_img_docker),
        alt   = title,
        style = "margin:0px auto;display:block"
    )

    return(network_screenshot)
}


get_em_view <- function(title, network = NULL) {
    # Purpose:
    #     Export a screenshot of the current network on a host
    # Parameters:
    #     title: Title of the network
    #     host_dir: working directory of the host machine
    #     network: Cytoscape suid of a network
    # Value:
    #     network_screenshot: screen shot of the EM
    
    if (!dir.exists(file.path(getwd(), "screenshots")))
        dir.create(file.path(getwd(), "screenshots"))
    
    fname <- gsub("[^a-zA-Z0-9]", "_", title)
    fname <- paste0(fname, ".png")
    
    export_network_img <- file.path(getwd(), "screenshots", fname)

    if (file.exists(export_network_img))
        response <- file.remove(export_network_img)
    
    suid     <- if (is.null(network)) NULL else as.numeric(network)
    response <- RCy3::exportImage(export_network_img,
                                  network = suid,
                                  type = "png")
    
    network_screenshot <- htmltools::img(
        src   = knitr::image_uri(export_network_img),
        alt   = title,
        style = "margin:0px auto;display:block"
    )
    
    return(network_screenshot)
}

# ====  TESTS  =================================================================
# Function tests to be added here...


# [END]
