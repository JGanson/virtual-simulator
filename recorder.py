import sounddevice as sd
from scipy.io.wavfile import write


fs = 44100
sec = 5
device = 0

i = 0
while True:
    answer = input("Record? ('Y'es/'N'o/'C'onfigure): ").upper()

    if answer == 'Y':
        print("Recording...")
        record = sd.rec(int(fs * sec), samplerate=fs, channels=2, device=device)
        sd.wait()

        write("recording" + str(i) + ".wav", fs, record)
        print("Recording saved!")

        i += 1
    elif answer == 'N':
        break
    elif answer == 'C':
        print("1. Recorder Configuration")
        print("2. Recording Time Configuration")
        print("3. Quit")

        config_option = input("Configuration Option: ").upper()

        if config_option == "1":
            print(sd.query_devices())
    else:
        print("Your answer is not valid. Please try again.")