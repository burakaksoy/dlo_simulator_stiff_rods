#!/bin/bash
sleep 1s;

# gnome-terminal --tab --title="ROSCORE" --command "bash -c \"source ~/.bashrc; killall gzclient && killall gzserver; roscore; exec bash\"";
# sleep 1s;

gnome-terminal --tab --title="Rviz" --command "bash -c \"source ~/.bashrc; roslaunch dlo_simulator_stiff_rods rviz.launch; exec bash\"";
sleep 4s;

gnome-terminal --tab --title="Simulator" --command "bash -c \"source ~/.bashrc; roslaunch dlo_simulator_stiff_rods dlo_simulator.launch; exec bash\"";
sleep 1s;

gnome-terminal --tab --title="Test GUI" --command "bash -c \"source ~/.bashrc; rosrun dlo_simulator_stiff_rods test_gui.py; exec bash\"";
sleep 1s;

# Use it to publish spacenav_twist msg
gnome-terminal --tab --title="RQT EZ Pub" --command "bash -c \"source ~/.bashrc; rosrun rqt_ez_publisher rqt_ez_publisher; exec bash\"";
sleep 1s;




