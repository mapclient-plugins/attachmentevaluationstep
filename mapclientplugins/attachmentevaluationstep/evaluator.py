"""
MAP Client, a program to generate detailed musculoskeletal models for OpenSim.
    Copyright (C) 2012  University of Auckland
    
This file is part of MAP Client. (http://launchpad.net/mapclient)

    MAP Client is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MAP Client is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MAP Client.  If not, see <http://www.gnu.org/licenses/>..
"""

"""
Class for evaluating attachment site coordinates on a fieldwork model
"""

import os
import numpy as np
from gias.common import simplemesh_tools
from gias.musculoskeletal.bonemodels.modelcore import BoneAttachmentRegions
from fieldwork.field import geometric_field

SELF_DIRECTORY = os.path.split(__file__)[0]

ATTACHMENT_FILES = {
                    'femur - left': os.path.join(SELF_DIRECTORY, 'data/FEMUR_LEFT_ATTACHMENT_ATLAS24x24.xml'),
                    'tibia - left': os.path.join(SELF_DIRECTORY, 'data/TIBIA_ATTACHMENT_UPD_ATLAS24x24.xml'),
                    # 'fibula - left': os.path.join(SELF_DIRECTORY, 'data/'),
                    # 'patella - left': os.path.join(SELF_DIRECTORY, 'data/'),
                    # 'hemipelvis - left': os.path.join(SELF_DIRECTORY, 'data/'),
                    }
VALID_MODELS = set(ATTACHMENT_FILES.keys())

class Evaluator(object):
    """ Class for evaluating the coordinates of attachment sites on a model
    """

    def __init__(self, **kwargs):
        self.model = None
        self._modelsm = None
        self.attachment_coordinates = {}
        self._model_name = None
        self._regions = None
        self._vertices = None
        self._faces = None

    @property
    def model_name(self):
        return self._model_name

    @model_name.setter
    def model_name(self, name):
        if name not in VALID_MODELS:
            raise ValueError('invalid model name')

        self._model_name = name
        self._load_attachment_file()

    def _load_attachment_file(self):
        self._regions = BoneAttachmentRegions()
        self._regions.from_xml(ATTACHMENT_FILES[self.model_name])
    
    def evaluate_attachment_coordinates(self):
        disc = [int(x) for x in self._regions.atlas_disc.split('x')]
        self._vertices, self._faces = self.model.triangulate(disc, merge=True)

        self.attachment_coordinates = {}
        for ri, r in enumerate(self._regions.regions):
            self.attachment_coordinates['_'.join([r.name, r.end, r.number])] = self._vertices[r.vertices,:].tolist()

    def make_labelled_mesh(self):

        mesh = simplemesh_tools.simpleMesh(v=self._vertices, f=self._faces)

        vert_labels = np.zeros(len(self._vertices), dtype=int)
        # face_labels = np.zeros(len(self._faces), dtype=int)
        code = {}
        for ri, r in enumerate(self._regions.regions):
            vert_labels[r.vertices] = ri+1
            # face_labels[r.faces] = ri+1
            code[ri+1] = (r.name, r.end, r.number)

        return mesh, vert_labels, face_labels, code