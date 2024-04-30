import math
from particle import Particle

class ForceModel:
    def apply(self, particle):
        # 应用力到颗粒上，需要在子类中具体实现
        pass

# 首先定义 Force 类的子类，分别对应重力、弹性力、阻尼力和粘附力

class GravityForce(ForceModel):
    def __init__(self, g):
        self.g = g

    def apply(self, particle):
        # 计算重力
        F = particle.mass * self.g
        return F

class ElasticForce(ForceModel):
    def __init__(self):
        pass

    def apply(self, particle):
        # 计算弹性力
        # 总为斥力，仅和坐标相关
        if particle.distance_contact >= 0:
            return 0
        else:
            F = particle.direction_contact * particle.spring_coeff * abs(particle.distance_contact)
            print(particle.direction_contact, particle.spring_coeff, abs(particle.distance_contact))
            return F

class DampingForce(ForceModel):
    def __init__(self):
        pass

    def apply(self, particle):
        # 计算阻尼力
        # 总与接近速度方向相反，颗粒速度-壁面速度
        if particle.distance_contact >= 0:
            return 0
        else:
            F = -1 * particle.damping_coeff * particle.velocity_contact
            return F

class AdhesiveForce(ForceModel):
    def __init__(self):
        pass
    
    def apply(self, particle):
        # 计算粘附力
        # 总为引力，仅和坐标有关
        if particle.distance_contact > Particle.adhesive_outer_cutoff:
            F_ad = 0
        else:
            if particle.distance_contact > Particle.adhesive_inner_cutoff:
                F_ad = particle.Hamaker_contact * particle.radius_eff / (6*particle.distance_contact**2) * \
                    (particle.rough / (particle.rough + particle.radius_eff) + 1/((1+particle.rough/particle.distance_contact)**2))
            else:
                F_ad = 4 * math.pi * particle.surface_energy_contact * particle.radius_eff * (particle.rough / \
                    (particle.rough + particle.radius_eff) + 1/((1+particle.rough/Particle.adhesive_inner_cutoff)**2))
        # 粘附力方向保持向下
        return -1 * particle.direction_contact * F_ad