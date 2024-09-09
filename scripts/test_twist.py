#!/usr/bin/env python3
import rospy
from geometry_msgs.msg import Twist


def publish_velocity():
    rospy.init_node('custom_velocity_publisher', anonymous=True)

    # Publisher for Twist messages
    # pub = rospy.Publisher('/custom_static_particles_cmd_vel_0', Twist, queue_size=1)
    pub = rospy.Publisher('/cmd_vel_particle_0', Twist, queue_size=1)
    pub2 = rospy.Publisher('/cmd_vel_particle_39', Twist, queue_size=1)
    
    rate = rospy.Rate(100)  # 10 Hz

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
        total_elapsed = (current_time - start_time).to_sec()

        # Velocity update logic
        if total_elapsed < zero_velocity_duration:
            linear_velocity = 0.0
        elif total_elapsed < zero_velocity_duration + ramp_up_duration:
            linear_velocity = (total_elapsed - zero_velocity_duration) / ramp_up_duration * max_velocity
        elif total_elapsed < ramp_down_start:
            linear_velocity = max_velocity
        else:
            # Ramp down logic
            ramp_down_end = ramp_down_start + ramp_down_duration
            if total_elapsed < ramp_down_end:
                linear_velocity = max_velocity * (1.0 - (total_elapsed - ramp_down_start) / ramp_down_duration)
            else:
                linear_velocity = 0.0

        # Prepare and publish Twist message
        twist_msg = Twist()
        twist_msg.linear.x = linear_velocity
        # Set angular velocities if needed, e.g., twist_msg.angular.z = angular_velocity

        pub.publish(twist_msg)
        pub2.publish(twist_msg)
        last_time = current_time
        rate.sleep()

if __name__ == '__main__':
    try:
        publish_velocity()
    except rospy.ROSInterruptException:
        pass
