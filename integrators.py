from copy import deepcopy
import numpy as np

class Integrator:
    def __init__(self):
        pass

    def integrate(self, particle, forces, dt):
        # 根据选择的积分方法更新颗粒的状态
        pass

class RungeKuttaIntegrator(Integrator):
    pass
#     def integrate(self, particle, particles, grounds, forces, collision_forces, dt):
#         # 复制颗粒群以创建虚拟颗粒群
#         virtual_particles = deepcopy(particles)

#         # 计算所有颗粒的k1速度和位置
#         k1_pos = np.array([[p.particleID ,p.position] for p in virtual_particles])
#         k1_vel = np.array([[p.particleID ,p.velocity] for p in virtual_particles])
#         k1_acc = np.array([[p.particleID ,p.current_force / p.mass] for p in virtual_particles])
        
#         k2_pos = k1_pos + k1_vel*dt/2
#         k2_vel = k1_vel + k1_acc*dt/2
        
#         # 更新虚拟颗粒群状态以计算k2
#         for i, vp in enumerate(virtual_particles):
#             vp.velocity = k1_vel[i] + k1_acc[i] * dt / 2
#             vp.position = k1_pos[i] + k1_vel[i] * dt / 2
#         self.update_virtual_particles(virtual_particles, grounds, forces, collision_forces)
#         k2_vel = np.array([vp.velocity for vp in virtual_particles])
#         k2_pos = np.array([vp.position for vp in virtual_particles])
#         k2_acc = np.array([vp.current_force / vp.mass for vp in virtual_particles])

#         # 类似地计算k3和k4

#         # 使用k1, k2, k3, k4更新实际颗粒群的状态
#         for i, particle in enumerate(particles):
#             particle.velocity += (k1_vel[i] + 2 * k2_vel[i] + 2 * k3_vel[i] + k4_vel[i]) / 6
#             particle.position += (k1_pos[i] + 2 * k2_pos[i] + 2 * k3_pos[i] + k4_pos[i]) / 6

class EulerIntegrator(Integrator):
    def integrate(self, particle, particles, grounds, forces, collision_forces, dt):
        # 实现 Euler 积分方法
        acc = sum(particle.current_force) / particle.mass
        particle.velocity += acc * dt
        particle.position += particle.velocity * dt