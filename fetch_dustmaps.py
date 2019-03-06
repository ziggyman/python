from dustmaps.config import config
config['data_dir'] = '/Volumes/obiwan/azuri/data/dustmaps'

import dustmaps.sfd
dustmaps.sfd.fetch()

import dustmaps.planck
dustmaps.planck.fetch()

import dustmaps.bayestar
dustmaps.bayestar.fetch()

import dustmaps.iphas
dustmaps.iphas.fetch()

import dustmaps.marshall
dustmaps.marshall.fetch()

import dustmaps.chen2014
dustmaps.chen2014.fetch()

import dustmaps.lenz2017
dustmaps.lenz2017.fetch()

import dustmaps.pg2010
dustmaps.pg2010.fetch()
