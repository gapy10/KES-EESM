import matplotlib.pyplot as plt
import numpy as np

def plt_kpt(uind):
    tok=[]
    for i in range(len(uind)):
        tok.append((i)*5)

    plt.figure(figsize=(8, 5))
    plt.plot(tok, uind, marker='o')
    plt.xlabel('Vzbujalni tok [A]')
    plt.ylabel('Efektivna inducirana napetost [V]')
    plt.title('Efektivna inducriana napetost v odvisnosti od vzbujalnega toka')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plt_uind(uind):
    kot=[]
    for i in range(len(uind)):
        kot.append((i+1)*5)

    plt.figure(figsize=(8, 5))
    plt.plot(kot, uind, marker='o')
    plt.xlabel('Kolesni kot [stopinje]')
    plt.ylabel('Inducirana napetost [V]')
    plt.title('Inducirana napetost v odvisnosti od kolesnega kota')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    torques=torques = np.array(uind)
    angles = np.array(kot)
    N = len(torques)
    torques -= np.mean(torques)

    # FFT
    torque_fft = np.fft.fft(torques)
    magnitudes = 2 * np.abs(torque_fft[:N // 2]) / N

    # Frequency axis
    d_angle = angles[1] - angles[0]
    frequencies_deg = np.fft.fftfreq(N, d=d_angle)
    f_mech = 133.33  # Hz
    frequencies_Hz = frequencies_deg[:N // 2] * f_mech * 360

    # Extract harmonic components (multiples of 400 Hz)
    harmonic_base = 400  # electrical fundamental frequency
    num_harmonics = 10
    target_freqs = [harmonic_base * i for i in range(1, num_harmonics + 1)]

    # Find nearest bins and collect magnitudes
    harmonic_amplitudes = []
    actual_freqs = []

    for target in target_freqs:
        idx = np.argmin(np.abs(frequencies_Hz - target))
        actual_freqs.append(frequencies_Hz[idx])
        harmonic_amplitudes.append(magnitudes[idx])

    print("U1: ", harmonic_amplitudes[0])
    print("U3: ", harmonic_amplitudes[2])
    # Plot
    plt.figure(figsize=(10, 5))
    plt.bar(range(1, num_harmonics + 1), harmonic_amplitudes, tick_label=[f"{int(f)} Hz" for f in target_freqs])
    plt.xlabel("Harmonska komponenta (Hz)")
    plt.ylabel("Amplituda [V]")
    plt.title("Harmonske komponente inducirane napetosti")
    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.show()

def plt_M(M):
    kot=[]
    for i in range(len(M)):
        kot.append(i*5)

    plt.figure(figsize=(8, 5))
    plt.plot(kot, M, marker='o')
    plt.xlabel('Kolesni kot [stopinje]')
    #plt.ylabel('Navor [Nm]')
    plt.ylabel('Navor [Nm]')
    plt.title('Navor v odvisnosti od kolesnega kota')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    torques=torques = np.array(M)
    angles = np.array(kot)
    N = len(torques)
    torques -= np.mean(torques)

    # FFT
    torque_fft = np.fft.fft(torques)
    magnitudes = 2 * np.abs(torque_fft[:N // 2]) / N

    # Frequency axis
    d_angle = angles[1] - angles[0]
    frequencies_deg = np.fft.fftfreq(N, d=d_angle)
    f_mech = 133.33  # Hz
    frequencies_Hz = frequencies_deg[:N // 2] * f_mech * 360

    # Extract harmonic components (multiples of 400 Hz)
    harmonic_base = 400  # electrical fundamental frequency
    num_harmonics = 10
    target_freqs = [harmonic_base * i for i in range(1, num_harmonics + 1)]

    # Find nearest bins and collect magnitudes
    harmonic_amplitudes = []
    actual_freqs = []

    for target in target_freqs:
        idx = np.argmin(np.abs(frequencies_Hz - target))
        actual_freqs.append(frequencies_Hz[idx])
        harmonic_amplitudes.append(magnitudes[idx])

    print("M1: ",harmonic_amplitudes[0])
    print("M3: ",harmonic_amplitudes[2])
    # Plot
    plt.figure(figsize=(10, 5))
    plt.bar(range(1, num_harmonics + 1), harmonic_amplitudes, tick_label=[f"{int(f)} Hz" for f in target_freqs])
    plt.xlabel("Harmonska komponenta (Hz)")
    plt.ylabel("Amplituda [Nm]")
    plt.title("Harmonske komponente navora na rotor")
    plt.grid(True, axis='y')
    plt.tight_layout()
    plt.show()

