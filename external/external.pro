TEMPLATE = subdirs
SUBDIRS  = quadprog++    \
           quadprog++-v2

system("python --help") {

	system("python obtain_eigen.py")
	system("python obtain_ceres.py")
	system("python obtain_nlopt.py")
	system("python obtain_ipopt.py")

} else {
	warning("Python has not been found, you will need to compile external dependancies on your own.")
}