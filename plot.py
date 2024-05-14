import matplotlib.pyplot as plt


def task_b():
    with open('data/B/gaus_end.txt', 'r') as file:
        gaus_file = file.readlines()
        gaus_norm_history = []
        for line in gaus_file:
            gaus_norm_history.append(float(line))

    with open('data/B/jaco_end.txt', 'r') as file:
        jaco_file = file.readlines()
        jaco_norm_history = []
        for line in jaco_file:
            jaco_norm_history.append(float(line))

    plt.plot(jaco_norm_history, label="Jacobi")
    plt.plot(gaus_norm_history, label="Gauss Seidel")
    plt.yscale('log')
    plt.xlabel('Iterations')
    plt.ylabel('Value')
    plt.title('Residual error history')
    plt.grid(True)
    plt.legend()
    plt.savefig('figure0.png') 
    plt.show()

def task_c():
    with open('data/B/gaus_inf.txt', 'r') as file:
        gaus_file = file.readlines()
        gaus_norm_history = []
        for line in gaus_file:
            gaus_norm_history.append(float(line))

    with open('data/B/jaco_inf.txt', 'r') as file:
        jaco_file = file.readlines()
        jaco_norm_history = []
        for line in jaco_file:
            jaco_norm_history.append(float(line))
    plt.figure()
    plt.plot(jaco_norm_history, label="Jacobi")
    plt.plot(gaus_norm_history, label="Gauss Seidel")
    plt.yscale('log')
    plt.xlabel('Iterations')
    plt.ylabel('Value')
    plt.title('Residual error history')
    plt.grid(True)
    plt.legend() 
    plt.savefig('figure1.png')
    plt.show()
   

def task_e():
    with open('data/E/time_direct.txt', 'r') as file:
        direct_file = file.readlines()
        direct_times = []
        for line in direct_file:
            direct_times.append(float(line))

    with open('data/E/time_jacobs.txt', 'r') as file:
        jaco_file = file.readlines()
        jaco_times = []
        for line in jaco_file:
            jaco_times.append(float(line))
    with open('data/E/time_gaus.txt', 'r') as file:
        gaus_file = file.readlines()
        gaus_times = []
        for line in gaus_file:
            gaus_times.append(float(line))

    sizes = [100, 500, 1000, 2000, 3000, 4000, 5000, 6000]

    plt.figure()
    plt.plot(sizes, direct_times, label="LU")
    plt.plot(sizes, jaco_times, label="Jacobi")
    plt.plot(sizes, gaus_times, label="Gauss Seidel")
    plt.xlabel('Sizes')
    plt.ylabel('Time in seconds')
    plt.title('Execution time comparison')
    plt.grid(True)
    plt.legend() 
    plt.savefig('figure2.png')
    plt.show()
    
    pass

if __name__ == "__main__":
    task_b()
    task_c()
    task_e()
