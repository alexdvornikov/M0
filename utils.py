import numpy as np
import h5py
import yaml
import consts
from consts import detector

clock_interval = 0.1 # us

def hit_to_3d(geometry, hits, last_trigger):

    hits_pos = [[], [], []]
    hits_index = []
    hits_charge = []
    hits_text = []

    for ih, hit in enumerate(hits):

        io_group = hit['iogroup']
        io_channel = hit['iochannel']
        chip = hit['chipid']
        channel = hit['channelid']

        module_id = (io_group - 1)//4

        x_offset = geometry.tpc_offsets[module_id][0]*10
        y_offset = geometry.tpc_offsets[module_id][1]*10
        z_offset = geometry.tpc_offsets[module_id][2]*10

        io_group = io_group - (io_group - 1)//4*4
        x,y = geometry.geometry[(io_group, io_channel, chip, channel)]
        
        hits_pos[0].append(x + x_offset)
        hits_pos[1].append(y + y_offset)
        hits_pos[2].append(geometry.get_z_coordinate(io_group, io_channel, clock_interval*(hit['ts'] - last_trigger)) + z_offset)

    return hits_pos

def get_TPC_bounds():
    bounds = []
    for ix in range(detector.TPC_BORDERS.shape[0]):
        bounds.append(detector.TPC_BORDERS[ix]*10)
    return np.array( bounds )

def get_track_ends(track, my_geometry):
    track_start = np.array([track['start'][0] + my_geometry.tpc_offsets[0][0]*10,
                            track['start'][1] +my_geometry.tpc_offsets[0][1]*10,
                            track['start'][2] + my_geometry.tpc_offsets[0][2]*10])
    track_end = np.array([track['end'][0] + my_geometry.tpc_offsets[0][0]*10,
                          track['end'][1] + my_geometry.tpc_offsets[0][1]*10,
                          track['end'][2] + my_geometry.tpc_offsets[0][2]*10])
    return track_start, track_end

def shift_track_to_anode(track_start, track_end, anode_z):

    p1_z = track_start[2]
    p2_z = track_end[2]

    if p1_z > 0: #If first TPC (anode at positive z)
        shift = anode_z - p1_z
        p1_z_shifted = anode_z
        p2_z_shifted = p2_z + shift

    if p1_z < 0: #If second TPC (anode at negative z)
        shift = anode_z - abs(p1_z) 
        p1_z_shifted = -anode_z
        p2_z_shifted = p2_z - shift

    return p1_z_shifted, p2_z_shifted

class DetectorGeometry():
    """Class describing the geometry of the Detector"""

    def __init__(self, detector_properties, pixel_layout):
        self.detector_properties = detector_properties
        self.pixel_layout = pixel_layout
        self.geometry = {}
        self.io_group_io_channel_to_tile = {}
        self.tile_positions = None
        self.tile_orientations = None
        self.tpc_offsets = None
        self.load_geometry()
        consts.load_properties(detector_properties,pixel_layout)

    @staticmethod
    def rotate_pixel(pixel_pos, tile_orientation):
        return pixel_pos[0]*tile_orientation[2], pixel_pos[1]*tile_orientation[1]

    def load_geometry(self):
        geometry_yaml = yaml.load(open(self.pixel_layout), Loader=yaml.FullLoader)

        pixel_pitch = geometry_yaml['pixel_pitch']
        chip_channel_to_position = geometry_yaml['chip_channel_to_position']
        self.tile_orientations = geometry_yaml['tile_orientations']
        self.tile_positions = geometry_yaml['tile_positions']
        tile_chip_to_io = geometry_yaml['tile_chip_to_io']
        xs = np.array(list(chip_channel_to_position.values()))[:, 0] * pixel_pitch
        ys = np.array(list(chip_channel_to_position.values()))[:, 1] * pixel_pitch
        x_size = max(xs)-min(xs)+pixel_pitch
        y_size = max(ys)-min(ys)+pixel_pitch

        tile_geometry = {}

        with open(self.detector_properties) as df:
            detprop = yaml.load(df, Loader=yaml.FullLoader)

        self.tpc_offsets = detprop['tpc_offsets']
        for tile in tile_chip_to_io:
            tile_orientation = self.tile_orientations[tile]
            tile_geometry[tile] = self.tile_positions[tile], self.tile_orientations[tile]

            for chip in tile_chip_to_io[tile]:
                io_group_io_channel = tile_chip_to_io[tile][chip]
                io_group = io_group_io_channel//1000
                io_channel = io_group_io_channel % 1000
                self.io_group_io_channel_to_tile[(io_group, io_channel)] = tile

            for chip_channel in chip_channel_to_position:
                chip = chip_channel // 1000
                channel = chip_channel % 1000

                try:
                    io_group_io_channel = tile_chip_to_io[tile][chip]
                except KeyError:
                    continue

                io_group = io_group_io_channel // 1000
                io_channel = io_group_io_channel % 1000
                x = chip_channel_to_position[chip_channel][0] * \
                    pixel_pitch - x_size / 2 + pixel_pitch / 2
                y = chip_channel_to_position[chip_channel][1] * \
                    pixel_pitch - y_size / 2 + pixel_pitch / 2

                x, y = self.rotate_pixel((x, y), tile_orientation)
                x += self.tile_positions[tile][2]
                y += self.tile_positions[tile][1]
                self.geometry[(io_group, io_channel, chip, channel)] = x, y
                
    def get_z_coordinate(self, io_group, io_channel, time):
        tile_id = self.get_tile_id(io_group, io_channel)

        z_anode = self.tile_positions[tile_id][0]
        drift_direction = self.tile_orientations[tile_id][0]
        return z_anode + time * detector.V_DRIFT * drift_direction * detector.cm/detector.mm

    def get_tile_id(self, io_group, io_channel):
        if (io_group, io_channel) in self.io_group_io_channel_to_tile:
            tile_id = self.io_group_io_channel_to_tile[io_group, io_channel]
        else:
            return np.nan

        return tile_id
