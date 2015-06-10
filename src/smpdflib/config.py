# -*- coding: utf-8 -*-
"""
Created on Tue May  5 15:03:41 2015

@author: zah
"""
from __future__ import print_function
import sys
import itertools
import glob
import fnmatch
from collections import Counter

import yaml

import smpdflib.lhaindex as lhaindex
import smpdflib.actions as actions
from smpdflib.core import PDF, make_observable

class ConfigError(ValueError): pass

class Config(object):
    def __init__(self, params):
        self._group_counter = iter("group_%d_" % i for i in itertools.count())
        self._grid_names = []


        self.actiongroups = self.parse_params(params)


    _group_len = None

    def parse_defaults(self, group):
        if not isinstance(group, dict):
            raise ConfigError("Defaults not understood: %s" % group)
        if 'prefix' in group:
            raise ConfigError("'prefix' not allowed in defaults")
        return group

    def parse_action_group(self,  group, defaults=None):
        if not isinstance(group, dict):
            raise ConfigError("Group not understood: %s" % group)

        if defaults is None:
            defaults = {}
        else:
            defaults = defaults.copy()

        d = {}
        if 'actions' in group:
            acts = group['actions']
        elif 'actions' in defaults:
            acts = defaults['actions']
        else:
            acts = {'all'}

        acts = self.parse_actions(acts)

        if 'observables' in group:
            observables = self.parse_observables(group['observables'])
        elif 'observables' in defaults:
            observables = self.parse_observables(defaults['observables'])
        else:
            if any(actions.requires_result(act) for act in acts):
                #plot_asq does not need observables.
                raise ConfigError("No observables found for action group.")
            observables = []

        if 'pdfsets' in group:
            pdfsets = self.parse_pdfsets(group['pdfsets'])
        elif 'pdfsets' in defaults:
            pdfsets = self.parse_pdfsets(defaults['pdfsets'])
        else:
            raise ConfigError("No pdfsets found for action group.")

        if 'base_pdf' in group or 'base_pdf' in defaults:
            if 'base_pdf' in group:
                base_pdf = self.parse_base_pdf(group['base_pdf'])
            elif 'base_pdf' in defaults:
                base_pdf =  self.parse_base_pdf(defaults['base_pdf'])
            if base_pdf not in pdfsets:
                raise ConfigError("Base pdf must be included in pdfsets: %s"
                                  % base_pdf)
            d['base_pdf'] = base_pdf

        if 'smpdf_spec' in group:
            smpdf_spec = self.parse_smpdf_spec(group['smpdf_spec'],
                                                     observables)
        elif 'smpdf_spec' in defaults:
            smpdf_spec = self.parse_smpdf_spec(defaults['smpdf_spec'],
                                                     observables)
        else:
            smpdf_spec = [{'Observables': observables}]



        groupcount = next(self._group_counter)
        if not 'prefix' in group:
            if self._group_len > 1:
                group['prefix'] = groupcount
            else:
                group['prefix'] = ''
        d['observables'] = observables
        d['pdfsets'] = pdfsets
        d['actions'] = acts
        d['smpdf_spec'] = smpdf_spec
        #Substitute the things we just parsed
        group.update(d)
        #Dump group on top of defaults
        defaults.update(group)
        #Now it has everything
        final = defaults
        #Finally make the final checks, that require everything to be parsed.
        self.check_actiongroup(final)
        return final

    def check_actiongroup(self, final):
        for action in final['actions']:
            actionfunc = actions.ACTION_DICT[action]
            if hasattr(actionfunc, 'appends_to'):
                final[actionfunc.appends_to] = []
            if hasattr(actionfunc, 'checks'):
                try:
                    for check_func in actionfunc.checks:
                        check_func(action, final, self)
                except actions.ActionError as e:
                    raise ConfigError(e.message)


    def parse_pdfsets(self, pdfs):
        pdfsets =  []
        for pdf in pdfs:
            if isinstance(pdf, dict):
                try:
                    names = pdf['name']
                except KeyError:
                    raise ConfigError("Unrecognized "
                                      "format for pdfsets: %s" % pdfs)
            elif isinstance(pdf, str):
                names = pdf
            else:
                raise ConfigError("Unrecognized format for pdfsets: %s" % pdfs)

            newsets = [PDF(name) for name in
                        (lhaindex.expand_local_names(names) +
                        fnmatch.filter(self._grid_names, names))]
            if not newsets:
                raise ConfigError("pdfset is empty for specification '%s'. "
                "Is it in the LHAPDF path?" % names)
            remote_sets = {PDF(name) for
                           name in lhaindex.expand_index_names(names)}
            diff = remote_sets - set(newsets)
            if diff:
                #TODO: Use logging.
                print("WARNING: The specification '%s' matches "
                      "some official PDF sets that are not locally installed. "
                      "They are:\n%s" % (names, '\n'.join(str(d) for
                                                          d in diff)),
                      file=sys.stderr)
            pdfsets += newsets
        return pdfsets

    def parse_base_pdf(self, base):
        if isinstance(base, dict):
            try:
                name = base['name']
            except KeyError:
                raise ConfigError("Unrecognized format for pdfsets: %s" % base)
        elif isinstance(base, str):
            name = base
        else:
            raise ConfigError("Unrecognized format for pdfsets: %s" % base)
        existing = lhaindex.expand_local_names(name)
        if not existing:
            raise ConfigError("No PDF set %s. Is it in the LHAPDF path?"%name)
        if len(existing) > 1:
            raise ConfigError("Only one base_pdf allowed. Matches were %s"
                              % existing)
        return PDF(name)



    def parse_actions(self, acts):
        try:
            result = actions.build_actions(acts)
        except ValueError as e:
            raise ConfigError("Could not parse actions '%s': %s" % (acts,
                                                                    e.message))
        return result

    def parse_observables(self, obslitst):
        observables = []
        allnames = []
        for obs in obslitst:
            if isinstance(obs, str):
                obsdict = {'name':obs}
            elif isinstance(obs, dict):
                obsdict = obs.copy()
            else:
                raise ConfigError("Observable format not understood: %s" % obs)
            names = glob.glob(obsdict.pop('name'))
            if not names:
                raise ConfigError("No observables found for %s" %
                                  names)
            allnames += names
        for name in allnames:
            try:
                obsobj = make_observable(name, **obsdict)
            except ValueError as e:
                raise ConfigError("Could not parse the observable %s: %s"
                                       % (name, e.message))
            except TypeError as e:
                    raise ConfigError("Incorrect arguments passed to process "
                                      "observable %s: %s" % (name, e.message))
            observables.append(obsobj)
        s = set(observables)
        if len(s) != len(observables):
            #TODO: proper logging
            c = Counter(observables)
            dups = [obs.name for obs in observables if c[obs]>1]
            obs = c.keys()
            print("WARNING: Duplicate observables: %s" % dups, file=sys.stderr)

        return observables

    def parse_smpdf_spec(self, smpdf_spec, observables):
        if not isinstance(smpdf_spec, list):
            raise ConfigError("smpdf_spec must be a list of:\n"
                              "- prefix: {obslist}")
        d = {}
        all_obs = set()
        for item in smpdf_spec:
            if len(item) != 1:
                raise ConfigError("smpdf_spec must be a list of:\n"
                  "- prefix: {obslist}")
            key, obslist = item.items()[0]
            obshere = self.parse_observables(obslist)
            #Do not duplicate objects, to make use of caches and so on
            d[key] = [obs for obs in observables if obs in obshere]
            if any(obs not in observables for obs in obshere):
                raise ConfigError("Observable %s in smpdf_spec must be "
                                  "also in obsevables")
            s = set(obshere)
            common = all_obs & s
            if common:
                raise ConfigError("Duplicate observables in different "
                "smpdf_specs: %s" % [obs.name for obs in common])
            all_obs |= s


        return d



    def parse_params(self, params):
        if not 'actiongroups' in params:
            actiongroups = [params]
            defaults = {}
        else:
            actiongroups = params.pop('actiongroups')
            self._group_len = len(actiongroups)
            defaults = self.parse_defaults(params)
        actiongroups = [self.parse_action_group(group, defaults)
                        for group in actiongroups]
        return actiongroups

    @classmethod
    def from_yaml(cls, stream):
        try:
            return cls(yaml.load(stream))
        except yaml.error.MarkedYAMLError as e:
            raise ConfigError("Failed to parse yaml file: %s" % e)