from pynwb import NWBHDF5IO
import matplotlib.pyplot as plt
import numpy as np

def read_nwb(file_path):
    io = NWBHDF5IO(file_path, 'r')
    nwbfile = io.read()
    nwbfile._io = io

    print("📂 Información del archivo NWB")
    print("Descripción de sesión:", nwbfile.session_description)
    print("Institución:", nwbfile.institution)
    print("Experimento:", nwbfile.experiment_description)

    electrodes = list(nwbfile.icephys_electrodes)
    print("Electrodos intracelulares:", electrodes)

    acquisitions = list(nwbfile.acquisition)
    stimuli = list(nwbfile.stimulus)
    print("Adquisiciones disponibles:", acquisitions)
    print("Estímulos disponibles:", stimuli)

    raw_sweeps = [k for k in acquisitions if "raw_current" in k]
    print(f"Total de sweeps crudos: {len(raw_sweeps)}")

    return nwbfile, raw_sweeps


def close_nwb(nwbfile):
    io = getattr(nwbfile, "_io", None)
    if io is not None:
        io.close()

def plot_sweep(nwbfile, sweep_name):
    if sweep_name not in nwbfile.acquisition:
        raise ValueError(f"Sweep {sweep_name} no encontrado")
    
    sweep = nwbfile.acquisition[sweep_name]
    
    data = np.asarray(sweep.data[:])
    
    if hasattr(sweep, 'timestamps') and sweep.timestamps is not None:
        t = np.array(sweep.timestamps)
    else:
        t = np.arange(len(data)) / sweep.rate
    
    plt.figure(figsize=(8, 4))
    plt.plot(t, data)
    plt.title(sweep_name)
    plt.xlabel("Tiempo (s)")
    plt.ylabel(f"Corriente [{sweep.unit}]")
    plt.show()
