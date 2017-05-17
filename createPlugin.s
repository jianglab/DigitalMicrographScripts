
//The following simple script gives an example of the command syntax for plug-in creation. There are 7 arguments in the command.

//1st argument: When this scipt is executed in DM, it looks at the given path for the given script file (*.s), e.g. "C:\\Documents and Settings\\vdbad\\Desktop\\s2stigmator_CM200_DM1_plugin\\Astigmatism_MagDistortion20_RY_GUI.s"
//Note the "\\" in the directory

//2nd argument: If this script (e.g. Astigmatism_MagDistortion20_RY_GUI.s) is found, it will generate a plug-in with the given name "s2stigmator" (e.g. s2stigmator.gt1)
//Note that different DM versions generate the plug-in in different folders: 
//	        For DM version 1.x, C:\Program Files\Gatan\DigitalMicrograph\Plugins
//	        For DM version 2.x or 3.x, C:\Users\YourUsername\AppData\Local\Gatan\Plugins

//3rd argument: The first number 0 in the example definds the "package level" can can be 0, 1, 2 or 3. 
//The value of this parameter sets the plug-in file extension as *.gt1, *.gt2, *.gt3 or *.gtk, respectively.
//This controls the order in which plug-ins are loaded when DM is launched. This is important where one plug-in makes calls on functions in a second plug-in.

//4th, 5th and 6th argument:
//              If this script is found, it will install it as a menu command with the given name (4th argument), 
//	        in the given menu (5th argument)
//	        and the optional submenu (6th argument)

//7th argument: The second number 0 can be 0 or 1. This value dicates whether the script wil be installed as a menu command or a library, respectively.

//In general, the command is:
//              AddScriptFileToPackage(your-script-and-directory , plugin-name-in-generation , package-level , plugin-name-in-menu, menu-name , submenu-name , 0)

//When installing the plug-in to other computers, copy the s2stigmator plug-in into the DM plug-in-folder:
//              For DM version 1.x, C:\Program Files\Gatan\DigitalMicrograph\PlugIns\
//              For DM version 2.x or 3.x, C:\Program Files\Gatan\PlugIns\

//For more information, please visit https://www.felmi-zfe.at/cms/wp-content/uploads/2016/07/how_to_organize_dm_scripts.pdf

AddScriptFileToPackage( "C:\\Documents and Settings\\vdbad\\Desktop\\s2stigmator_CM200_DM1_plugin\\Astigmatism_MagDistortion20_RY_GUI.s" , "s2stigmator" , 0 , "s2stigmator" , "s2stigmator" , "" , 0)