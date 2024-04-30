import json
import pandas as pd
from particle import Particle, Ground
from forces import GravityForce, ElasticForce, DampingForce, AdhesiveForce
from integrators import *
from simulation import Simulation
import matplotlib.pyplot as plt

def load_config(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

def main(file_path):
    # Load configuration
    config = load_config(file_path)

    # Load material properties
    df_materials = pd.read_excel(config["material_properties_file"])
    df_materials.set_index('Properties', inplace=True)
    
    # 设置类变量
    Particle.adhesive_inner_cutoff = config.get('adhesive_inner_cutoff')
    Particle.adhesive_outer_cutoff = config.get('adhesive_outer_cutoff')  

    # Initialize particles
    particles = [Particle(particleID=p["particleID"], radius=p["radius"],
                          material_name=p["material_name"],
                          initial_velocity=p["initial_velocity"],
                          initial_height=p["initial_height"], df_materials=df_materials, rough=p["roughness"])
                 for p in config["particles"]]

    # Initialize simulation
    sim = Simulation(dt=config["time_step"], total_time=config["total_simulation_time"])

    # Add particles to simulation
    for particle in particles:
        sim.add_particle(particle)

    # Create and add ground
    ground_config = config["ground"]
    ground = Ground(groundID= ground_config["groundID"], material_name=ground_config["material_name"], 
                    position=ground_config["position"], df_materials=df_materials, rough=ground_config["roughness"])
    sim.add_ground(ground)

    # Initialize force models
    gravity_force = GravityForce(g=config["gravity"])
    elastic_force = ElasticForce()
    damping_force = DampingForce()
    adhesive_force = AdhesiveForce()

    # Add forces to simulation
    sim.add_force(gravity_force)
    sim.add_force(elastic_force,is_collision_force=True)
    sim.add_force(damping_force,is_collision_force=True)
    sim.add_force(adhesive_force,is_collision_force=True)

    # Set integrator
    integrator_type = config["integrator"]
    if integrator_type == "RungeKutta":
        integrator = RungeKuttaIntegrator()
    elif integrator_type == "Euler":
        integrator = EulerIntegrator()
    else:
        raise ValueError("Unknown integrator type: " + integrator_type)
    sim.set_integrator(integrator)

    # Run simulation
    particleResults = sim.run()
    
    return particleResults

if __name__ == "__main__":
    file_path = './particle_config.json'
    particleResults = main(file_path)
    
    # 假设 sim.particles 是包含两个 Particle 实例的列表
    # 假设每个 Particle 的 history 属性是一个包含 [time, position, force] 的列表

    # 提取第一个颗粒的时间和位置数据
    time1, position1, force1, velocity1 = zip(*[(record[0], record[1], record[2], record[3]) for record in particleResults[0].history])
    
    # 确定最大长度
    force1 = list(force1)
    max_length = max(len(item) for item in force1)
    forces1 = [row + [np.nan] * (max_length - len(row)) for row in force1]
    forces1 = np.array(forces1)
    
    # # 提取第二个颗粒的时间和位置数据
    # time2, position2 = zip(*[(record[0], record[1]) for record in self.particles[1].history])

    incident_velocity = velocity1[0]
    rebound_velocity = velocity1[-1]
    COR = rebound_velocity / abs(incident_velocity)
    print(f"入射速度：{incident_velocity:.2f} m/s ;  回弹速度：{rebound_velocity:.2f} m/s")
    print(f"恢复系数 : {COR:.2f}")

    # 绘制散点图
    plt.scatter(time1, position1, label='Particle 1 Position')
    # plt.scatter(time2, position2, label='Particle 2 Position')

    # 设置图表标题和坐标轴标签
    plt.title('Particle Position Over Time')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.legend()
    
    # 绘制散点图
    plt.figure(2)
    plt.scatter(time1, velocity1, label='Particle 1 velocity')

    # 设置图表标题和坐标轴标签
    plt.title('Particle Velocity Over Time')
    plt.xlabel('Time')
    plt.ylabel('Velocity')
    plt.legend()

    # 绘制散点图
    plt.figure(3)
    forcesName = ['gravity_force', 'elastic_force', 'damping_force', 'adhesive_force']
    for i in range(forces1.shape[1]):
        plt.scatter(time1, forces1[:,i], label=forcesName[i])  # 绘制每一列，添加图例标签
    # 设置图表标题和坐标轴标签
    plt.title('Particle force Over Time')
    plt.xlabel('Time')
    plt.ylabel('force')
    plt.legend()
    
    # 显示图表
    plt.show()
        