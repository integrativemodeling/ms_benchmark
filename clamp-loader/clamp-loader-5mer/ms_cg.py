#!/usr/bin/python

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
                               IMP.algebra.Vector3D(200, 200, 200))


# this function creates the molecular hierarchies for the various involved
# proteins
def create_representation():
    m = IMP.Model()
    all = IMP.atom.Hierarchy.setup_particle(IMP.Particle(m))
    all.set_name("the universe")
    # create a protein, represented as a set of connected balls of appropriate
    # radii and number, chose by the resolution parameter and the number of
    # amino acids.

    def create_protein(name, ds):
        h = IMP.atom.create_protein(m, name, resolution, ds)
        leaves = IMP.atom.get_leaves(h)
        # for convenience, have one molecular hierarchy containing all
        # molecules
        all.add_child(h)
        r = IMP.atom.create_connectivity_restraint([IMP.atom.Selection(c)
                                                    for c in h.get_children()],
                                                   k)
        if r:
            m.add_restraint(r)
            # only allow the particles to penetrate or separate by 1 angstrom
            m.set_maximum_score(r, k)
    create_protein("ProteinG1", 431)
    create_protein("ProteinG2", 431)
    create_protein("ProteinG3", 431)
    create_protein("ProteinD", 343)
    create_protein("ProteinDP", 334)
#    create_protein("ProteinCHI", 147)
#    create_protein("ProteinPSI", 137)

 #   create_protein("ProteinD", 110)
 #   create_protein("ProteinE", 125)
    return (m, all)


# create the needed restraints and add them to the model
def create_restraints(m, all):
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
        iG = r.add_type([rps[0], rps[1], rps[2]])
        iD = r.add_type([rps[3]])
        iDP = r.add_type([rps[4]])
        #iCHI = r.add_type([rps[5]])
        #iPSI = r.add_type([rps[6]])
#        iD = r.add_type([rps[5]])

        n1 = r.add_composite([iG, iG, iG, iD, iDP])
        n2 = r.add_composite([iG, iG, iG], n1)
        n3 = r.add_composite([iG, iG], n2)

        n4 = r.add_composite([iG, iG, iG, iDP], n1)
        n5 = r.add_composite([iG, iG, iDP], n4)
        n6 = r.add_composite([iG, iDP], n5)
        n7 = r.add_composite([iD, iDP], n1)

# Data from other sources: binary interacrions (x-linkis, affinity purificationP
#        n3 = r.add_composite([iA, iA], n2)
#        n5 = r.add_composite([iA, iA, iC, iD, iE], n4)
#        n7 = r.add_composite([iA, iA, iA, iD, iE], n4)
#        n8 = r.add_composite([iA, iA, iA, iC], n4)
#        n12 = r.add_composite([iB, iC], n11)
#        n6 = r.add_composite([iA, iA, iD, iE], n5)
#        n6a = r.add_composite([iA, iA, iD, iE], n7)
#        n9 = r.add_composite([iA, iA, iC], n8)
#        n10 = r.add_composite([iA, iC], n9)
        m.add_restraint(r)
        m.set_maximum_score(r, k)

    def add_distance_restraint(s0, s1):
        r = IMP.atom.create_distance_restraint(s0, s1, 0, k)
        m.add_restraint(r)
        # only allow the particles to separate by one angstrom
        m.set_maximum_score(r, k)
    # evr=IMP.atom.create_excluded_volume_restraint([all])

    evr = IMP.atom.create_excluded_volume_restraint([all])

    #evr= IMP.core.ExcludedVolumeRestraint(all,1)
    m.add_restraint(evr)
    # a Selection allows for natural specification of what the restraints act
    # on
    S = IMP.atom.Selection
    sG1 = S(hierarchy=all, molecule="ProteinG1")
    sG2 = S(hierarchy=all, molecule="ProteinG2")
#    sA3=S(hierarchy=all, molecule="ProteinA3")
    sG3 = S(hierarchy=all, molecule="ProteinG3")
    sD = S(hierarchy=all, molecule="ProteinD")
    sDP = S(hierarchy=all, molecule="ProteinDP")
#    sCHI=S(hierarchy=all, molecule="ProteinCHI")
#    sPSI=S(hierarchy=all, molecule="ProteinPSI")

    add_connectivity_restraint([sG1, sG2, sG3, sD, sDP])  # , sPSI, sCHI])
    # for l in IMP.atom.get_leaves(all):
    #    r= IMP.example.ExampleRestraint(l, k)
    #    m.add_restraint(r)
    #    m.set_maximum_score(k)
#    print m.evaluate(False)

# find acceptable conformations of the model


def get_conformations(m):
    sampler = IMP.core.MCCGSampler(m)
    sampler.set_bounding_box(bb)
    # magic numbers, experiment with them and make them large enough for
    # things to work
    sampler.set_number_of_conjugate_gradient_steps(250)
    sampler.set_number_of_monte_carlo_steps(50)
    sampler.set_number_of_attempts(10000)
    # We don't care to see the output from the sampler
    sampler.set_log_level(IMP.SILENT)

    # return the IMP.ConfigurationSet storing all the found configurations that
    # meet the various restraint maximum scores.
    cs = sampler.get_sample()
    return cs

# cluster the conformations and write them to a file


def analyze_conformations(cs, all, gs):
    # we want to cluster the configurations to make them easier to understand
    # in the case, the clustering is pretty meaningless
    embed = IMP.statistics.ConfigurationSetXYZEmbedding(cs,
                                                        IMP.container.ListSingletonContainer(IMP.atom.get_leaves(all)), True)
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
(m, all) = create_representation()
IMP.atom.show_molecular_hierarchy(all)
create_restraints(m, all)


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

cs = get_conformations(m)
# print cs

# Report solutions
print "found", cs.get_number_of_configurations(), "solutions"

for i in range(0, cs.get_number_of_configurations()):
    cs.load_configuration(i)
    # print the configuration
    print "solution number: ", i, "scored :", m.evaluate(False)


ListScores = []
for i in range(0, cs.get_number_of_configurations()):
    cs.load_configuration(i)
    # print the configuration
    print "solution number: ", i, "scored :", m.evaluate(False)
    ListScores.append(m.evaluate(False))
    print ListScores

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

analyze_conformations(cs, all, gs)
# print m.evaluate(False)
