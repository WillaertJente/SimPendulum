# SimPendulum
Simulations for pendulum test 

The master folder contains 3 folders (BK, Implicit and Forward)

BK:      
The BK folder contains experimental data to track (BK/DF0) and personal information on the subject. 

Forward:    
The forward folder contains the script to run the pendulum simulations forward.                     
This script works for the torque model with and without SRS.                         
This script also works for the muscle model with and without SRS, but not yet for the model with feedback (work in progress)                    

Implicit:   
Implicit-torque:   
This script tracks pendulum kinematics based on a torque model with damping, SRS and feedback                                  
Implicit-muscle:    
This folder contains a seperate folder 'Experimental data', which is the experimental data that should be used for the muscle model.                           
Pendulum_spiermodel_v1 is a simulation without SRS and works                          
Pendulum_spiermodel_v2_SRS is still work in progress. The goal is to add SRS.                             
