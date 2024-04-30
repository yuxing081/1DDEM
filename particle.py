import math
import numpy as np
from fun import Function

class Particle:
    # 将表格列名映射到类的属性名
    property_map = {
        'Young modulus (Pa)': 'young_modulus',
        'Density (kg/m3)': 'density',
        'Poisson ratio': 'poisson_ratio',
        'Surface energy (J/m2)': 'surface_energy',
        'Hamaker (J)': 'hamaker'
        # 可以根据需要添加更多映射
    }
    
    # 初始的类变量定义，可以设为默认值
    adhesive_inner_cutoff = None
    adhesive_outer_cutoff = None


    def __init__(self, particleID, radius, material_name, df_materials, rough, initial_velocity, spring_coeff = 1000, initial_height=0, COE_contact = 0.8):
        self.particleID = particleID
        self.radius = radius
        self.load_material_properties(material_name, df_materials)
        self.rough = rough
        self.velocity = initial_velocity
        self.spring_coeff = spring_coeff
        self.position = initial_height
        self.COE_contact = COE_contact
                
        self.mass = self.calculate_mass()
        self.history = []  # 添加用于记录受力历史的字典
        self.current_force = []  # 添加用于存储当前时间步受力的属性

    def load_material_properties(self, material_name, df_materials):
        material_properties = df_materials.loc[material_name]
        for column, attribute in self.property_map.items():
            setattr(self, attribute, material_properties[column])
        # 根据需要可以添加更多属性的设置

    def calculate_mass(self):
        volume = (4.0 / 3.0) * math.pi * self.radius**3
        return self.density * volume
    
    def check_collision_with_ground(self, ground):
        distance = math.sqrt((self.position - ground.position)**2)
        return distance < self.radius + self.adhesive_outer_cutoff

    def check_collision_with_particle(self, other):
        distance = math.sqrt((self.position - other.position)**2)
        return distance < self.radius + other.radius + self.adhesive_outer_cutoff
    
        # 新增方法：记录
    def record(self, time):
        current_record = [time, self.position, self.current_force, self.velocity]
        self.history.append(current_record)

    # 新增方法：重置当前受力
    def reset_force(self):
        self.current_force = []

    # 修改方法：应用力并累加到当前受力
    def apply_force(self, force):
        self.current_force.append(force)
     
    def calculate_contact_properties(self, other_particle, dt):
        # 以表面能为例，实际计算应根据您的模型和需求进行
        self.surface_energy_contact = (self.surface_energy * other_particle.surface_energy) ** 0.5
        self.Hamaker_contact = (self.hamaker * other_particle.hamaker) ** 0.5
        # 计算有效半径和质量
        self.radius_eff = Function.Eff_phy_quantity(self.radius, other_particle.radius)
        self.mass_contact = Function.Eff_phy_quantity(self.mass, other_particle.mass)
        self.young_contact = 1/((1-self.poisson_ratio**2)/self.young_modulus+(1-other_particle.poisson_ratio**2)/other_particle.young_modulus)
        # 计算相对位置和速度
        self.distance_contact = abs(self.position - other_particle.position) - (self.radius + other_particle.radius)
        self.velocity_contact = self.velocity - other_particle.velocity
        self.direction_contact = (self.position - other_particle.position) / abs(self.position - other_particle.position)
        # 根据相对速度获取恢复系数,从而计算阻尼系数
        self.contact_time = 2.868*(((self.mass_contact**2)/((self.young_contact**2)*self.radius_eff*abs(self.velocity_contact))))**(0.2)
        if dt > self.contact_time or dt < self.contact_time/50:
            dt = self.contact_time/50
        self.spring_coeff = self.mass_contact/(self.contact_time**2) * ((np.log(self.COE_contact))**2+math.pi**2)
        if self.COE_contact == 0:
            self.damping_coeff = 2*(self.spring_coeff*self.mass_contact)**0.5
        else:
            self.damping_coeff = 2*(self.spring_coeff*self.mass_contact)**0.5 * abs(np.log(self.COE_contact)) / \
                (math.pi**2+(np.log(self.COE_contact))**2)**0.5
        return dt
        
    def calculate_wall_contact_properties(self, other_particle, dt):
        # 以表面能为例，实际计算应根据您的模型和需求进行
        self.surface_energy_contact = (self.surface_energy * other_particle.surface_energy) ** 0.5
        self.Hamaker_contact = (self.hamaker * other_particle.hamaker) ** 0.5
        # 计算有效半径和质量
        self.radius_eff = self.radius
        self.mass_contact = self.mass
        self.young_contact = 1/((1-self.poisson_ratio**2)/self.young_modulus+(1-other_particle.poisson_ratio**2)/other_particle.young_modulus)
        # 计算相对位置和速度
        self.distance_contact = abs(self.position - other_particle.position) - self.radius
        self.velocity_contact = self.velocity - other_particle.velocity
        self.direction_contact = (self.position - other_particle.position) / abs(self.position - other_particle.position)
        # 根据相对速度获取恢复系数,从而计算阻尼系数
        # 注意相对速度为负值时会报错        
        self.contact_time = 2.868*(((self.mass_contact**2)/((self.young_contact**2)*self.radius_eff*abs(self.velocity_contact))))**(0.2)
        if dt > self.contact_time or dt < self.contact_time/50:
            dt = self.contact_time/50
        self.spring_coeff = self.mass_contact/(self.contact_time**2) * ((np.log(self.COE_contact))**2+math.pi**2)
        if self.COE_contact == 0:
            self.damping_coeff = 2*(self.spring_coeff*self.mass_contact)**0.5
        else:
            self.damping_coeff = 2*(self.spring_coeff*self.mass_contact)**0.5 * abs(np.log(self.COE_contact)) / \
                (math.pi**2+(np.log(self.COE_contact))**2)**0.5
        return dt

    def calculate_total_force(self, other_particles, grounds, forces, collision_forces, dt):
        # 重置当前受力
        self.reset_force()

        # 应用普通力（如重力、风力等）
        for force in forces:
            self.apply_force(force.apply(self))
            

        # 检查并处理与地面的碰撞
        for ground in grounds:
            if self.check_collision_with_ground(ground):
                dt = self.calculate_wall_contact_properties(ground,dt)
                for collision_force in collision_forces:
                    self.apply_force(collision_force.apply(self))

        # 检查并处理与其他颗粒的碰撞
        for other_particle in other_particles:
            if self.particleID != other_particle.particleID and self.check_collision_with_particle(other_particle):
                dt = self.calculate_contact_properties(other_particle,dt)
                for collision_force in collision_forces:
                    self.apply_force(collision_force.apply(self))
        return dt
            
class Ground(Particle):
    def __init__(self, groundID, material_name, df_materials, position, rough):
        # 正确传递参数给父类的构造方法
        super().__init__(particleID=-1, radius=float('inf'), material_name=material_name, df_materials=df_materials, initial_height=position, rough=rough,initial_velocity=0)
        self.groundID = groundID
