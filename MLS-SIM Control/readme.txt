Code
\Supplementary Software\MLS-SIM SystemControl

Launcher.vi
This program is the launcher for the whole system control program. This program could run directly when the camera and NI cards are connected to the computer directly (see 'wiring guide' below).
Prerequisites:
LabVIEW, 64-bit, version 2021 or higher.
8GB RAM

MLS-SIM.lvproj
The project file of this LabVIEW project. This project is constructed with the actor framework.

Main_Control\
This folder contains the actor controlling the main user interface and its related functions.

Camera\
This folder contains the actor controlling the camera and its related functions.

NI\
This folder contains the actor controlling NI cards and its related functions.

Utilities\
This folder contains functions used by all the actors.

MainControl FP Params.ini
This file is an inner parameter file of the program and cannot be directly modified. It can only be modified through the program itself. This file contains the front panel control parameters in the main control UI. These parameters define all the imaging progress such as frame size, exposure time, frame interval time, color mode, etc.

Camera Inner Parameters.ini
This file is an inner parameter file of the program and cannot be directly modified. It can only be modified through the program itself. This file contains the innate parameters of the connected camera. These parameters contain the full field of view size and the output trigger port index.

Cam FP Params.ini
This file is an inner parameter file of the program and cannot be directly modified. It can only be modified through the program itself. This file contains the front panel control parameters in the image display UI. These parameters define all the display progress such as binning number, brightness, contrast, etc.

NI Inner Parameters.ini
This file is an inner parameter file of the program and cannot be directly modified. It can only be modified through the program itself. This file contains the front panel control parameters in the phase monitor UI.

Wiring guide:
The trigger output port 2 of the Hamamatsu camera is connected to /Dev0/PFI0, which triggers the following tasks:
1. FastClock task: Counter output. Trigger is /Dev0/PFI0. Terminal is /Dev0/PFI2. Clock is inner clock.
2. MainTrig task: Counter output. Terminal is /Dev0/PFI1. Clock is /Dev0/PFI0.
3. AO1 laser task: Analog output. Trigger is /Dev0/PFI1. Terminal links to the AOTF control. Clock is /Dev0/PFI2.
4. AO2 zy task: Analog output. Trigger is /Dev0/PFI1. Terminal links to the piezo and the scanning mirror control. Clock is /Dev0/PFI0.
5. AI phase task: Analog input. Clock is /Dev0/PFI0.
