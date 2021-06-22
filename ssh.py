import paramiko




def main():
    """Description of main() - what does this function do?  Does it run a 
		program?  Does it execute test code?"""

    host = "brccluster.cos.ncsu.edu"
    port = 22
    username = "yyang55"
    password = "1xuanYang@nc"
    commands = "sbatch -p bigmem -n 64 tcoffeeCommand.sh\n"#CreateCommand()


    print("commands are uploading")
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, port, username, password)
    stdin, stdout, stderr = ssh.exec_command(commands)
    lines = stdout.readlines()
    for line in lines:
        print(line[:-1])



if __name__ == '__main__':
    main()
