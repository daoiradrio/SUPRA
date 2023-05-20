from time import sleep



print("[", end="", flush=True)
for i in range(50):
    print("=", end="", flush=True)
    sleep(0.1)
print("]")
