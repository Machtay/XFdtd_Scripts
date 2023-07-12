# XFdtd_Scripts
Some basic scripts for setting up geometries and simulations in Remcom's XFdtd.

The directory Example_Scripts will hold some examples I've made or modified from things given to me by XF. Anyone with access to the repository is welcome to add some example scripts there that might be helpful for people to look at.

The GENETIS_Scripts directory contains xmacro scripts used by GENETIS. There's some room for cleanup there that's actually currently happening on one of the GENETIS repositories. For practice purposes, the script simulation_PEC.xmacro is probably the best place to start, since it is comprised of everything you need to set up a simulation. Other than changing a path to the data, it should be good to work out of the box.

XF.conf is my config file. I haven't done anything fancy to it, it's just the default I seem to have. Sometimes people have a config file that causes bugs or doesn't look at clean, so feel free to copy this one into your XFdtd config directory.

xfsolver_array.sh is a bash script for submitting a batch job to OSC. OSC uses slurm for managing job submissions. This should also be good to run out of the box, provided you have an XF project set up with a simulation queued. There are instructions on how to run at the top: you'll need to replaced the variables as appropraite (NPOP should be the simulation number you want to run and XFProj=<path/to/XF/project>).

To get started, run the following commands in the GENETIS_Examples directory (replacing <path/to/your/clone> appropriately). These two commands will replace the path to my copy of the repo with the path to your copy:

vim -c ':%s/\/fs\/project\/PAS0654\/Machtay_XF_Tests\/Scripts\//<path/to/your/clone>' + -c ':wq!' output.xmacro

vim -c ':%s/\/fs\/project\/PAS0654\/Machtay_XF_Tests\/Scripts\//<path/to/your/clone>' + -c ':wq!' simulation_PEC.xmacro



