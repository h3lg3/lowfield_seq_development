"""Functions to write and read sequence definitions."""

from pypulseq.Sequence.sequence import Sequence

custom_seq_definitons =  ('fov', 'slice_thickness', 'Name', 
                          'Flipangle', 'number_of_readouts', 'k_space_encoding1', 'Nr', 'k_space_encoding2', 
                          'slices', 'average', 'phase', 'contrast', 
                          'repetition', 'set', 'segment', 'N_interleaves', 
                          'delta', 'sampling_scheme', 'TE', 'TR', 'proj_mode', 'train_duration',
                          'n_total_trains', 'tr_delay', 'channel_order', 'etl', 'ro_bandwidth', 'ro_oversampling')

def write_seq_definitions(
    seq: Sequence,
    fov: tuple,
    slice_thickness: float,
    name: str,
    alpha: float,
    Nx: int,
    Ny: int = 1,
    Nr: int = 1,
    Nz: int = 1,
    N_slices: int = 1,
    average: float = 1,
    phase: float = 1,
    contrast: float = 1,
    repetition: float = 1,
    set: float = 1,  # noqa: A002
    segment: float = 1,
    N_interleaves: float = 1,
    delta: float = 0,
    sampling_scheme: str = 'cartesian',
    TE: float = 0,
    TR: float = 0,
    proj_mode: bool = False,
    train_duration: float = 0,
    n_total_trains: int = 0,
    tr_delay: float = 0,
    channel_order: tuple = (),
    etl: int = 0,
    ro_bandwidth: float = 0,
    ro_oversampling: int = 1,
) -> None:
    """Write sequence definitions into the sequence object."""
    if sampling_scheme not in ['radial', 'cartesian', 'spiral']:
        raise TypeError('Unknown sampling scheme')

    seq.set_definition('FOV', fov)
    seq.set_definition('slice_thickness', slice_thickness)
    seq.set_definition('name', name)
    seq.set_definition('Flipangle', alpha)
    seq.set_definition('number_of_readouts', int(Nx))
    seq.set_definition('ro_oversampling', int(ro_oversampling))

    if sampling_scheme == 'radial' and Nr != 1:
        seq.set_definition('number_of_spokes', int(Nr))
    elif sampling_scheme == 'spiral' and N_interleaves != 1:
        seq.set_definition('N_interleaves', N_interleaves)
    else:
        seq.set_definition('k_space_encoding1', int(Ny))

    if N_slices != 1:
        seq.set_definition('slices', N_slices)
    if average != 1:
        seq.set_definition('average', average)
    if phase != 1:
        seq.set_definition('phase', phase)
    if contrast != 1:
        seq.set_definition('contrast', contrast)
    if repetition != 1:
        seq.set_definition('repetition', repetition)
    if set != 1:
        seq.set_definition('set', set)
    if segment != 1:
        seq.set_definition('segment', segment)
    if Nz != 1:
        seq.set_definition('k_space_encoding2', int(Nz))
    if delta != 0:
        seq.set_definition('delta', delta)

    seq.set_definition('sampling_scheme', sampling_scheme)
    if TE != 0:
        seq.set_definition('TE', TE)
    if TR != 0:
        seq.set_definition('TR', TR)

    if proj_mode:
        seq.set_definition('proj_mode', str(proj_mode))

    if train_duration != 0:
        seq.set_definition('train_duration', train_duration)
    if n_total_trains != 0:
        seq.set_definition('n_total_trains', n_total_trains)
    if tr_delay != 0:
        seq.set_definition('tr_delay', tr_delay)
    if channel_order != ():
        seq.set_definition('channel_order', channel_order)
    if etl != 0:
        seq.set_definition('etl', etl)        
    if ro_bandwidth != 0:
        seq.set_definition('ro_bandwidth', ro_bandwidth)           


def read_definitions(seq):
    """Read sequence definitions from a sequence object."""
    # default values
    Ny = 1
    Nr = 1
    Nz = 1
    N_slices = 1
    average = 1
    phase = 1
    contrast = 1
    repetition = 1
    set = 1  # noqa: A001
    segment = 1
    N_interleaves = 1
    TE = 0
    TR = 0
    alpha = 0
    proj_mode = False

    if 'FOV' in seq.dict_definitions:
        fov_x = seq.dict_definitions['FOV'][0]
        fov_y = seq.dict_definitions['FOV'][1]
        slice_thickness = seq.dict_definitions['FOV'][2]
    else:
        raise TypeError('FOV not given')

    if 'sampling_scheme' in seq.dict_definitions:
        Sampling_scheme = seq.dict_definitions['sampling_scheme'][0]

    if 'number_of_readouts' in seq.dict_definitions:  # kx
        Nx = seq.dict_definitions['number_of_readouts'][0]
    else:
        raise TypeError('number_of_readouts not given')

    if 'k_space_encoding1' in seq.dict_definitions:
        Ny = seq.dict_definitions['k_space_encoding1'][0]

    if 'number_of_spokes' in seq.dict_definitions:
        Nr = seq.dict_definitions['number_of_spokes'][0]

    if 'slices' in seq.dict_definitions:
        N_slices = seq.dict_definitions['slices'][0]

    if 'k_space_encoding2' in seq.dict_definitions:
        Nz = seq.dict_definitions['k_space_encoding2'][0]

    if 'average' in seq.dict_definitions:
        average = seq.dict_definitions['average'][0]

    if 'TE' in seq.dict_definitions:
        TE = seq.dict_definitions['TE'][0]

    if 'TR' in seq.dict_definitions:
        TR = seq.dict_definitions['TR'][0]

    if 'phase' in seq.dict_definitions:
        phase = seq.dict_definitions['phase'][0]

    if 'contrast' in seq.dict_definitions:
        contrast = seq.dict_definitions['contrast'][0]

    if 'repetition' in seq.dict_definitions:
        repetition = seq.dict_definitions['repetition'][0]

    if 'set' in seq.dict_definitions:
        set = seq.dict_definitions['set'][0]  # noqa: A001

    if 'segment' in seq.dict_definitions:
        segment = seq.dict_definitions['segment'][0]

    if 'N_interleaves' in seq.dict_definitions:
        N_interleaves = seq.dict_definitions['N_interleaves'][0]

    if 'Flipangle' in seq.dict_definitions:
        alpha = seq.dict_definitions['Flipangle'][0]

    if 'proj_mode' in seq.dict_definitions:
        proj_mode = bool(seq.dict_definitions['proj_mode'][0])

    dico = {
        'fov_x': fov_x,
        'fov_y': fov_y,
        'TE': TE,
        'TR': TR,
        'slice_thickness': slice_thickness,
        'sampling_scheme': Sampling_scheme,
        'Flip_angle': alpha,
        'number_of_readouts': round(Nx),
        'k_space_encoding1': round(Ny),
        'number_of_spokes': round(Nr),
        'slices': round(N_slices),
        'average': round(average),
        'phase': round(phase),
        'contrast': round(contrast),
        'repetition': round(repetition),
        'set': round(set),
        'segment': round(segment),
        'N_interleaves': round(N_interleaves),
        'k_space_encoding2': round(Nz),
        'proj_mode': proj_mode,
    }

    if dico.get('k_space_encoding1') == 1:
        dico['k_space_encoding1'] = dico['number_of_spokes']
    if dico.get('number_of_spokes') == 1:
        dico['number_of_spokes'] = dico['k_space_encoding1']

    if dico.get('k_space_encoding1') == 1 and dico.get('number_of_spokes') == 1:
        dico['k_space_encoding1'] = dico['number_of_readouts']

    return dico
