# -*- coding: utf-8 -*-
"""
Created on Tue May  5 15:03:41 2015

@author: zah
"""
import itertools
import glob

import yaml

import lhaindex
import actions
from smpdflib.core import PDF, Observable

class ConfigError(ValueError): pass

class Config(object):
    def __init__(self, actiongroups):
        self.actiongroups = actiongroups

    _group_counter = iter("group_%d_" % i for i in itertools.count())

    @classmethod
    def parse_defaults(cls, group):
        if not isinstance(group, dict):
            raise ConfigError("Defaults not understood: %s" % group)
        if 'observables' in group:
            group['observables'] = cls.parse_observables(group['observables'])
        if 'pdfsets' in group:
            group['pdfsets'] = cls.parse_pdfsets(group['pdfsets'])
        if 'actions' in group:
            acts = group['actions']
            group['actsions'] = actions.build_actions(acts)
        return group


    @classmethod
    def parse_action_group(cls,  group, defaults=None):
        if not isinstance(group, dict):
            raise ConfigError("Group not understood: %s" % group)

        if defaults is None:
            defaults = {}
        else:
            defaults = defaults.copy()

        d = {}
        if 'observables' in group:
            observables = cls.parse_observables(group['observables'])
        elif 'observables' in defaults:
            observables = defaults['observables']
        else:
            raise ConfigError("No observables found for action group.")

        if 'pdfsets' in group:
            pdfsets = cls.parse_pdfsets(group['pdfsets'])
        elif 'pdfsets' in defaults:
            pdfsets = defaults['pdfsets']
        else:
            raise ConfigError("No pdfsets found for action group.")

        if 'actions' in group:
            acts = group['actions']
        elif 'actions' in defaults:
            acts = defaults['actions']
        else:
            acts = {'all'}
        try:
            acts = actions.build_actions(acts)
        except ValueError as e:
            raise ConfigError(e.message)

        if 'base_pdf' in group:
            name = group['base_pdf']
            #TODO: Do dedicated method
            base = cls.parse_pdfsets([name])
            if len(base) > 1:
                raise ConfigError("Only one base allowed: %s" % name)
            base_pdf = base[0]
            if base_pdf not in pdfsets:
                raise ConfigError("Base pdf must be included in pdfsets: %s"
                                  % base_pdf)
            d['base_pdf'] = base_pdf

        groupcount = next(cls._group_counter)
        if not 'prefix' in group:
            group['prefix'] = groupcount
        d['observables'] = observables
        d['pdfsets'] = pdfsets
        d['actions'] = acts
        #Substitute the things we just parsed
        group.update(d)
        #Dump group on top of defaults
        defaults.update(group)
        #Now it has everything
        final = defaults
        return final


    @classmethod
    def parse_pdfsets(cls, pdfs):
        pdfsets =  []
        for pdf in pdfs:
            if isinstance(pdf, dict):
                names = pdf['name']
            elif isinstance(pdf, str):
                names = pdf
            else:
                raise ConfigError("Unrecognized format for pdfsets: %s" % pdfs)
            pdfsets += [PDF(name) for name in
                        lhaindex.expand_local_names(names)]
            if not pdfsets:
                raise ConfigError("pdfset is empty for specification '%s'. "
                "Is it in the LHAPDF path?"
                                  %names)
        return pdfsets

    @classmethod
    def parse_observables(cls, obslitst):
        observables = []
        for obs in obslitst:
             names = glob.glob(obs['name'])
             if not names:
                 raise ConfigError("No observables found for %s" % obs['name'])
             for name in names:
                 observables.append(Observable(name, obs['order']))
        return observables

    @classmethod
    def from_params(cls, params):
        if not 'actiongroups' in params:
            actiongroups = [params]
            defaults = {}
        else:
            actiongroups = params.pop('actiongroups')
            defaults = cls.parse_defaults(params)
        actiongroups = [cls.parse_action_group(group, defaults)
                        for group in actiongroups]
        return cls(actiongroups)

    @classmethod
    def from_yaml(cls, stream):
        try:
            return cls.from_params(yaml.load(stream))
        except yaml.parser.ParserError as e:
            raise ConfigError("Failed to parse yaml file: %s" % e)