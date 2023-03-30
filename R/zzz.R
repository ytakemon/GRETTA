.onAttach <- function(libname, pkgname) {
    hello <- r"{
    
    
    ██████╗ ██████╗ ███████╗████████╗ █████╗ 
    ██╔════╝ ██╔══██╗██╔════╝╚══██╔══╝██╔══██╗
    ██║  ███╗██████╔╝█████╗     ██║   ███████║
    ██║   ██║██╔══██╗██╔══╝     ██║   ██╔══██║
    ╚██████╔╝██║  ██║███████╗   ██║   ██║  ██║
    ╚═════╝ ╚═╝  ╚═╝╚══════╝   ╚═╝   ╚═╝  ╚═╝
    
    
    }"
  
    packageStartupMessage(
      cat(hello),
      "Welcome to GRETA! ", 
      "The version loaded is: ", utils::packageVersion("GRETA"), "\n",
      "The latest DepMap dataset accompanying this package is v22Q2. \n",
      "Please refer to our tutorial on GitHub for loading DepMap data and details: https://github.com/ytakemon/GRETA \n"
    )
}
