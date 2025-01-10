# wxsInit

Just a way to initialize variable for all scripts. You should edit paths to reference library and executables inside the functions. Yep, maybe not the simplest way but the easiest one to implement in my mind

- exec\_paths

    Send executable paths to the scripts this way. Directly edit the paths here

- data\_paths

    Send reference data paths to the scripts this way. Directly edit the paths here

- init\_conf

    Read init file and export a hash of variables. 

    Usage:

            my %wesconf = init_conf('project.init');
