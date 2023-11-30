.onAttach <- function(libname, pkgname) {
    hello <- r"{
      _______ .______       _______ .___________.___________.    ___      
     /  _____||   _  \     |   ____||           |           |   /   \     
    |  |  __  |  |_)  |    |  |__   `---|  |----`---|  |----`  /  ^  \    
    |  | |_ | |      /     |   __|      |  |        |  |      /  /_\  \   
    |  |__| | |  |\  \----.|  |____     |  |        |  |     /  _____  \  
     \______| | _| `._____||_______|    |__|        |__|    /__/     \__\ 
    
    }"
  
    packageStartupMessage(
      hello,
      "Welcome to GRETTA! ", 
      "The version loaded is: ", utils::packageVersion("GRETTA"), "\n",
      "The latest DepMap dataset accompanying this package is v23Q2. \n",
      "Please refer to our tutorial on GitHub for loading DepMap data and details: https://github.com/ytakemon/GRETTA \n"
    )
}
