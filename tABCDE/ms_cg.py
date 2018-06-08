#!/usr/bin/env python

from __future__ import print_function
import IMP
import IMP.atom
import IMP.container
import IMP.display
import IMP.statistics
import os
#import IMP.example

# the spring constant to use, it doesn't really matter
k = 100
# the target resolution for the representation, this is used to specify how detailed
# the representation used should be
resolution = 300
# the box to perform everything in
bb = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(0, 0, 0),
                               IMP.algebra.Vector3D(100, 100, 100))


# this function creates the molecular hierarchies for the various involved
# proteins
def create_representation():
    m = IMP.Model()
    rs = IMP.RestraintSet(m, "all restraints")
    all = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
    all.set_name("the universe")
    # create a protein, represented as a set of connected balls of appropriate
    # radii and number, chose by the resolution parameter and the number of
    # amino acids.

    def create_protein(name, ds):
        h = IMP.atom.create_protein(m, name, resolution, ds)
        leaves = IMP.atom.get_leaves(h)
        for l in leaves:
            IMP.core.XYZ(l).set_coordinates_are_optimized(True)
        # for convenience, have one molecular hierarchy containing all
        # molecules
        all.add_child(h)
        r = IMP.atom.create_connectivity_restraint([IMP.atom.Selection(c)
                                                    for c in h.get_children()],
                                                   k)
        if r:
            rs.add_restraint(r)
            # only allow the particles to penetrate or separate by 1 angstrom
            r.set_maximum_score(k)
    create_protein("ProteinA", 100)
    create_protein("ProteinB", 150)
    create_protein("ProteinC", 85)
    create_protein("ProteinD", 200)
    create_protein("ProteinE", 175)
 #   create_protein("ProteinD", 110)
 #   create_protein("ProteinE", 125)
    return (m, rs, all)


# create the needed restraints and add them to the restraint set
def create_restraints(m, rs, all):
    def add_connectivity_restraint(s):
        tr = IMP.core.TableRefiner()
        rps = []
        for sc in s:
            ps = sc.get_selected_particles()
            rps.append(ps[0])
            tr.add_particle(ps[0], ps)
        # duplicate the IMP.atom.create_connectivity_restraint functionality
        score = IMP.core.KClosePairsPairScore(
            IMP.core.HarmonicSphereDistancePairScore(0, 1),
            tr)
        # Work with IMP 2.1 or newer versions
        try:
            r = IMP.core.MSConnectivityRestraint(m, score)
        except NotImplementedError:
            r = IMP.core.MSConnectivityRestraint(score)
        iA = r.add_type([rps[0]])
        iB = r.add_type([rps[1]])
        iC = r.add_type([rps[2]])
        iD = r.add_type([rps[3]])
#        iD = r.add_type([rps[5]])
        iE = r.add_type([rps[4]])
        n1 = r.add_composite([iA, iB, iC, iD, iE])
        n2 = r.add_composite([iA, iB, iC, iD], n1)
        n3 = r.add_composite([iA, iB], n2)
        n4 = r.add_composite([iA, iB, iD], n2)
        n5 = r.add_composite([iA, iB, iE], n1)
        n6 = r.add_composite([iA, iB], n5)

# Data from other sources: binary interacrions (x-linkis, affinity
# purificationP
        n7 = r.add_composite([iB, iD], n1)
        n8 = r.add_composite([iA, iC], n1)

#        n3 = r.add_composite([iA, iA], n2)
#        n5 = r.add_composite([iA, iA, iC, iD, iE], n4)
#        n7 = r.add_composite([iA, iA, iA, iD, iE], n4)
#        n8 = r.add_composite([iA, iA, iA, iC], n4)
#        n12 = r.add_composite([iB, iC], n11)
#        n6 = r.add_composite([iA, iA, iD, iE], n5)
#        n6a = r.add_composite([iA, iA, iD, iE], n7)
#        n9 = r.add_composite([iA, iA, iC], n8)
#        n10 = r.add_composite([iA, iC], n9)
        rs.add_restraint(r)
        r.set_maximum_score(k)

    def add_distance_restraint(s0, s1):
        r = IMP.atom.create_distance_restraint(s0, s1, 0, k)
        rs.add_restraint(r)
        # only allow the particles to separate by one angstrom
        r.set_maximum_score(k)
    # evr=IMP.atom.create_excluded_volume_restraint([all])

    evr = IMP.atom.create_excluded_volume_restraint([all])

    #evr= IMP.core.ExcludedVolumeRestraint(all,1)
    rs.add_restraint(evr)
    # a Selection allows for natural specification of what the restraints act
    # on
    S = IMP.atom.Selection
    sA = S(hierarchy=all, molecule="ProteinA")
    sB = S(hierarchy=all, molecule="ProteinB")
#    sA3=S(hierarchy=all, molecule="ProteinA3")
    sC = S(hierarchy=all, molecule="ProteinC")
    sD = S(hierarchy=all, molecule="ProteinD")
    sE = S(hierarchy=all, molecule="ProteinE")
    add_connectivity_restraint([sA, sB, sC, sD, sE])
    # for l in IMP.atom.get_leaves(all):
    #    r= IMP.example.ExampleRestraint(l, k)
    #    rs.add_restraint(r)
    #    m.set_maximum_score(k)
#    print sf.evaluate(False)

# find acceptable conformations of the model


def get_conformations(m, sf):
    sampler = IMP.core.MCCGSampler(m)
    sampler.set_scoring_function(sf)
    sampler.set_bounding_box(bb)
    # magic numbers, experiment with them and make them large enough for
    # things to work
    sampler.set_number_of_conjugate_gradient_steps(200)
    sampler.set_number_of_monte_carlo_steps(50)
    sampler.set_number_of_attempts(10000)
    # We don't care to see the output from the sampler
    sampler.set_log_level(IMP.SILENT)

    # return the IMP.ConfigurationSet storing all the found configurations that
    # meet the various restraint maximum scores.
    cs = sampler.create_sample()
    return cs

# cluster the conformations and write them to a file


def analyze_conformations(m, cs, all, gs):
    # we want to cluster the configurations to make them easier to understand
    # in the case, the clustering is pretty meaningless
    embed = IMP.statistics.ConfigurationSetXYZEmbedding(cs,
                                                        IMP.container.ListSingletonContainer(m, IMP.atom.get_leaves(all)), True)
    cluster = IMP.statistics.create_lloyds_kmeans(embed, 10, 10000)
    # dump each cluster center to a file so it can be viewed.
    for i in range(cluster.get_number_of_clusters()):
        center = cluster.get_cluster_center(i)
        cs.load_configuration(i)
        w = IMP.display.PymolWriter("output/models/cluster.%d.pym" % i)
        for g in gs:
            w.add_geometry(g)


if not os.path.exists('output/models'):
    os.makedirs('output/models')

# now do the actual work
(m, rs, all) = create_representation()
IMP.atom.show_molecular_hierarchy(all)
create_restraints(m, rs, all)
sf = IMP.core.RestraintsScoringFunction([rs])


# in order to display the results, we need something that maps the particles onto
# geometric objets. The IMP.display.Geometry objects do this mapping.
# IMP.display.XYZRGeometry map an IMP.core.XYZR particle onto a sphere
gs = []
for i in range(all.get_number_of_children()):
    color = IMP.display.get_display_color(i)
    n = all.get_child(i)
    name = n.get_name()
    g = IMP.atom.HierarchyGeometry(n)
    g.set_color(color)
    gs.append(g)

cs = get_conformations(m, sf)
# print cs

# Report solutions
print("found", cs.get_number_of_configurations(), "solutions")

for i in range(0, cs.get_number_of_configurations()):
    cs.load_configuration(i)
    # print the configuration
    print("solution number: ", i, "scored :", sf.evaluate(False))


ListScores = []
for i in range(0, cs.get_number_of_configurations()):
    cs.load_configuration(i)
    # print the configuration
    print("solution number: ", i, "scored :", sf.evaluate(False))
    ListScores.append(sf.evaluate(False))
    print(ListScores)

f1 = open("output/scores.txt", "w")
f1.write("\n".join(map(lambda x: str(x), ListScores)))
f1.close()

# for each of the configuration, dump it to a file to view in pymol
for i in range(0, cs.get_number_of_configurations()):
    cs.load_configuration(i)
    w = IMP.display.PymolWriter("output/models/configuration.%d.pym" % i)
    for g in gs:
        w.add_geometry(g)
# Report solutions

analyze_conformations(m, cs, all, gs)
# print sf.evaluate(False)
