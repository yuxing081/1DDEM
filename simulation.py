import matplotlib.pyplot as plt

import numpy as np
import matplotlib.animation as animation

# 仿真类框架
class Simulation:
    def __init__(self, dt, total_time):
        self.particles = []
        self.grounds = []
        self.forces = []  # 用于存储体力，如重力
        self.collision_forces = []  # 用于存储碰撞力模型
        self.integrator = None
        self.dt = dt
        self.total_time = total_time

    def add_particle(self, particle):
        self.particles.append(particle)
        
    def add_ground(self, ground):
        self.grounds.append(ground)

    def add_force(self, force, is_collision_force=False):
        if is_collision_force:
            self.collision_forces.append(force)
        else:
            self.forces.append(force)

    def set_integrator(self, integrator):
        self.integrator = integrator

    def run(self):
        time = 0
        while time < self.total_time:
            # 计算并累加每个颗粒的受力
            for particle in self.particles:
                self.dt = particle.calculate_total_force(self.particles, self.grounds, self.forces, self.collision_forces, self.dt)
                # 记录每个颗粒的状态
                particle.record(time)

            # 使用积分器统一更新所有颗粒的状态
            for particle in self.particles:
                self.integrator.integrate(particle, self.particles, self.grounds, self.forces, self.collision_forces, self.dt)

            time += self.dt
        print("计算结束")
        
        return self.particles