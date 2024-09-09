#!/usr/bin/env python3
import rospy
from nav_msgs.msg import Odometry
from geometry_msgs.msg import Point, Quaternion
import tf.transformations as transformations
import numpy as np


def publish_odometry():
    rospy.init_node('custom_odom_publisher', anonymous=True)


    # pub = rospy.Publisher('/custom_static_particles_odom_0', Odometry, queue_size=1)
    pub = rospy.Publisher('/odom_particle_0', Odometry, queue_size=1)
    pub2 = rospy.Publisher('/odom_particle_39', Odometry, queue_size=1)
    rate = rospy.Rate(100)  # 10 Hz

    # position = Point(-1, 0, 2)
    # orientation = Quaternion(0, 0.7071081, 0, 0.7071055)
    
    position = Point(0.24, -0.08, 0.3)
    position2 = Point(-0.24, -0.08, 0.3)
    
    th_x, th_y,th_z = np.deg2rad([90, 0, -110])
    th_x2, th_y2,th_z2 = np.deg2rad([90, 0, -70])
    orientation = Quaternion(*transformations.quaternion_from_euler(th_x, th_y,th_z))
    orientation2 = Quaternion(*transformations.quaternion_from_euler(th_x2, th_y2,th_z2))

    # Time tracking
    start_time = rospy.Time.now()
    last_time = start_time

    # Phase durations in seconds
    max_velocity = 1
    zero_velocity_duration = 4
    ramp_up_duration = 0.1
    constant_velocity_duration = 0.4
    ramp_down_duration = 0.1
    ramp_down_start = zero_velocity_duration + ramp_up_duration + constant_velocity_duration

    while not rospy.is_shutdown():
        current_time = rospy.Time.now()
        elapsed_time = (current_time - last_time).to_sec()
        total_elapsed = (current_time - start_time).to_sec()

        # Velocity update logic
        if total_elapsed < zero_velocity_duration:
            velocity = 0.0
        elif total_elapsed < zero_velocity_duration + ramp_up_duration:
            velocity = (total_elapsed - zero_velocity_duration) / ramp_up_duration * max_velocity
        elif total_elapsed < ramp_down_start:
            velocity = max_velocity
        else:
            # Ramp down logic
            ramp_down_end = ramp_down_start + ramp_down_duration
            if total_elapsed < ramp_down_end:
                velocity = max_velocity * (1.0 - (total_elapsed - ramp_down_start) / ramp_down_duration)
            else:
                velocity = 0.0

        # Position integration
        position.x += velocity * elapsed_time
        position2.x += velocity * elapsed_time

        # Prepare and publish Odometry message
        odom_msg = Odometry()
        odom_msg.header.stamp = current_time
        odom_msg.header.frame_id = "odom"
        odom_msg.child_frame_id = "base_link"
        odom_msg.pose.pose.position = position
        odom_msg.pose.pose.orientation = orientation
        odom_msg.twist.twist.linear.x = velocity

        pub.publish(odom_msg)
        
        odom_msg2 = odom_msg
        odom_msg2.pose.pose.position = position2
        odom_msg2.pose.pose.orientation = orientation2
        odom_msg2.twist.twist.linear.x = velocity
        
        pub2.publish(odom_msg2)
        
        last_time = current_time
        rate.sleep()

if __name__ == '__main__':
    try:
        publish_odometry()
    except rospy.ROSInterruptException:
        pass
