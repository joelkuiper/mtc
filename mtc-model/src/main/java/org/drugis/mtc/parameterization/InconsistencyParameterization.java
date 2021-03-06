/*
 * This file is part of the GeMTC software for MTC model generation and
 * analysis. GeMTC is distributed from http://drugis.org/gemtc.
 * Copyright (C) 2009-2012 Gert van Valkenhoef.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.drugis.mtc.parameterization;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;

import org.drugis.mtc.graph.GraphUtil;
import org.drugis.mtc.graph.SpanningTreeIterable;
import org.drugis.mtc.model.Network;
import org.drugis.mtc.model.Study;
import org.drugis.mtc.model.Treatment;
import org.drugis.mtc.search.DepthFirstSearch;

import edu.uci.ics.jung.algorithms.transformation.FoldingTransformerFixed.FoldedEdge;
import edu.uci.ics.jung.graph.Hypergraph;
import edu.uci.ics.jung.graph.Tree;
import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.util.Pair;

/**
 * Implements parameterization of inconsistency models for network meta-analysis.
 * Method described in <a href="http://dx.doi.org/10.1007/s11222-011-9281-9">van Valkenhoef et al., Statistics and Computing (2011)</a>.
 */
public class InconsistencyParameterization extends ConsistencyParameterization {
	
	public static InconsistencyParameterization create(Network network) {
		Hypergraph<Treatment, Study> sGraph = NetworkModel.createStudyGraph(network);
		UndirectedGraph<Treatment, FoldedEdge<Treatment, Study>> cGraph = NetworkModel.createComparisonGraph(sGraph);
		Tree<Treatment, FoldedEdge<Treatment, Study>> tree = findSpanningTree(network.getStudies(), cGraph);
		Map<Partition, Set<List<Treatment>>> cycleClasses = getCycleClasses(cGraph, tree);
		Map<Study, Treatment> baselines = findStudyBaselines(network.getStudies(), cGraph, cycleClasses);
		InconsistencyParameterization pmtz = new InconsistencyParameterization(network, tree, cycleClasses, baselines);
		return pmtz;
	}
	
	/**
	 * Get the classes of fundamental cycles that share the same inconsistency
	 * factor (if any). The fundamental cycles are those cycles that are
	 * generated by the given spanning tree in the given study graph.
	 * (Note that the edge objects of the spanning tree need not equal those
	 * of the comparison graph, since the spanning tree features directed
	 * edges and the comparison graph is undirected. Thus any correspondence
	 * between the two needs to be established based on the vertices.)
	 * @param studyGraph The comparison graph.
	 * @param tree Spanning tree of the comparison graph.
	 * @return The cycle classes as identified by the reduced partitions.
	 */
	public static Map<Partition, Set<List<Treatment>>> getCycleClasses(UndirectedGraph<Treatment, FoldedEdge<Treatment, Study>> cGraph, Tree<Treatment, FoldedEdge<Treatment, Study>> tree) {
		Set<FoldedEdge<Treatment, Study>> nonTreeEdges = new HashSet<FoldedEdge<Treatment,Study>>(cGraph.getEdges());
		for (FoldedEdge<Treatment, Study> edge : tree.getEdges()) { // Go through vertices to remove tree edges from the set
			Pair<Treatment> vertices = new Pair<Treatment>(tree.getIncidentVertices(edge));
			nonTreeEdges.remove(cGraph.findEdge(vertices.getFirst(), vertices.getSecond()));
		}

		Map<Partition, Set<List<Treatment>>> cycleClasses = new HashMap<Partition, Set<List<Treatment>>>();

		// Each of the non-tree edges generates a fundamental cycle in the comparison graph.
		// We classify the fundamental cycles according to their reduced partitions.
		for (FoldedEdge<Treatment,Study> edge : nonTreeEdges) {
			Pair<Treatment> vertices = new Pair<Treatment>(cGraph.getIncidentVertices(edge));
			List<Treatment> cycle = GraphUtil.findPath(tree, vertices.getFirst(), vertices.getSecond());
			cycle.add(vertices.getFirst());
			cycle = standardizeCycle(cycle);
			
			List<Part> parts = new ArrayList<Part>(cycle.size() - 1);
			for (int i = 1; i < cycle.size(); ++i) {
				parts.add(new Part(cycle.get(i - 1), cycle.get(i), cGraph.findEdge(cycle.get(i - 1), cycle.get(i)).getFolded()));
			}
			Partition partition = new Partition(parts).reduce();
			
			if (!cycleClasses.containsKey(partition)) {
				cycleClasses.put(partition, new HashSet<List<Treatment>>());	
			}
			cycleClasses.get(partition).add(cycle);
		}
		
		return cycleClasses;
	}

	/**
	 * Determine the inconsistency degree given a classification of cycles into classes with equivalent reductions.
	 * @see getCycleClasses
	 */
	public static int getInconsistencyDegree(Map<Partition, Set<List<Treatment>>> cycleClasses) {
		int icd = 0;
		for (Partition p : cycleClasses.keySet()) {
			if (isInconsistencyCycle(p)) {
				++icd;
			}
		}
		return icd;
	}

	/**
	 * Determine whether the given partition represents an inconsistency cycle.
	 * @param p The partition.
	 * @return Whether the partition represents a potentially inconsistent cycle.
	 */
	public static boolean isInconsistencyCycle(Partition p) {
		return p.reduce().getParts().size() >= 3;
	}
	
	/**
	 * Find a baseline assignment that covers all comparisons.
	 * @param studies List of studies to find baselines for.
	 * @param cGraph Comparison graph.
	 */
	public static Map<Study, Treatment> findStudyBaselines(Collection<Study> studies, UndirectedGraph<Treatment, FoldedEdge<Treatment, Study>> cGraph) {
		return new DepthFirstSearch<Map<Study, Treatment>>().search(new InconsistencyBaselineSearchProblem(studies, cGraph));
	}
	
	/**
	 * Find a baseline assignment suitable for the given cycle classes.
	 * @param studies List of studies to find baselines for.
	 * @param cGraph Comparison graph.
	 * @param cycleClasses Cycle classes.
	 */
	public static Map<Study, Treatment> findStudyBaselines(Collection<Study> studies, UndirectedGraph<Treatment, FoldedEdge<Treatment, Study>> cGraph, Map<Partition, Set<List<Treatment>>> cycleClasses) {
		return new DepthFirstSearch<Map<Study, Treatment>>().search(new InconsistencyBaselineSearchProblem(studies, cGraph, cycleClasses));
	}
	
	public static Tree<Treatment, FoldedEdge<Treatment, Study>> findSpanningTree(Collection<Study> studies, UndirectedGraph<Treatment, FoldedEdge<Treatment, Study>> cGraph) {
		boolean hasCompleteBaseline = findStudyBaselines(studies, cGraph) != null;
		int nFunctional = cGraph.getEdgeCount() - cGraph.getVertexCount() + 1;
		int maxPossible = hasCompleteBaseline ? nFunctional : nFunctional - 1;
		
		int max = -1;
		Tree<Treatment, FoldedEdge<Treatment, Study>> best = null;
		Treatment root = CompareUtil.findLeast(cGraph.getVertices(), TreatmentComparator.INSTANCE);
		SpanningTreeIterable<Treatment, FoldedEdge<Treatment, Study>> spanningTreeIterable = new SpanningTreeIterable<Treatment, FoldedEdge<Treatment, Study>>(NetworkModel.toDirected(cGraph), root, TreatmentComparator.INSTANCE);
		for (Tree<Treatment, FoldedEdge<Treatment, Study>> tree : spanningTreeIterable) {
			Map<Partition, Set<List<Treatment>>> cycleClasses = getCycleClasses(cGraph, tree);
			int icd = getInconsistencyDegree(cycleClasses);
			if (icd > max && icd <= maxPossible) {
				if (hasCompleteBaseline || findStudyBaselines(studies, cGraph, cycleClasses) != null) {
					max = icd;
					best = tree;
				}
			}
			if (max == maxPossible) {
				return best;
			}
		}
		return best;
	}
	
	/**
	 * Standardize the given cycle by making the "least" treatment the first element,
	 * and by making the "least" neighbor of the first element the second element.
	 * @see TreatmentComparator
	 * @param cycle The cycle to standardize.
	 * @return Standardized version of the given cycle.
	 */
	public static List<Treatment> standardizeCycle(List<Treatment> cycle) {
		assertCycle(cycle);
		
		List<Treatment> std = rebaseCycle(cycle);
		
		// reverse the order if necessary
		if (TreatmentComparator.INSTANCE.compare(std.get(1), std.get(std.size() - 2)) > 0) {
			Collections.reverse(std);
		}
		
		return std;
	}

	private static List<Treatment> rebaseCycle(List<Treatment> cycle) {
		// find the least treatment
		Treatment least = CompareUtil.findLeast(cycle, TreatmentComparator.INSTANCE);
		
		// start the cycle from the least treatment
		List<Treatment> std = new ArrayList<Treatment>();
		std.addAll(cycle.subList(cycle.indexOf(least), cycle.size()));
		std.addAll(cycle.subList(1, cycle.indexOf(least) + 1));
		return std;
	}

	/**
	 * Check that the given list represents a cycle.
	 */
	private static void assertCycle(List<Treatment> cycle) {
		Set<Treatment> set = new HashSet<Treatment>(cycle);
		if (set.size() != cycle.size() - 1 || !cycle.get(0).equals(cycle.get(cycle.size() - 1))) {
			throw new IllegalStateException(cycle + " is not a cycle.");
		}
	}

	private final Map<Partition, Set<List<Treatment>>> d_cycleClasses;
	
	/**
	 * Construct an inconsistency parameterization with the given spanning tree and study baselines.
	 * @param network Network to parameterize.
	 * @param tree Spanning tree that defines the basic parameters.
	 * @param cycleClasses 
	 * @param baselines Map of studies to baseline treatments.
	 */
	public InconsistencyParameterization(Network network, Tree<Treatment, FoldedEdge<Treatment, Study>> tree, Map<Partition, Set<List<Treatment>>> cycleClasses, Map<Study, Treatment> baselines) {
		super(network, tree, baselines);
		d_cycleClasses = cycleClasses;
	}
	
	@Override
	public List<NetworkParameter> getParameters() {
		List<NetworkParameter> parameters = ConsistencyParameterization.getBasicParameters(d_tree);
		for (Entry<Partition, Set<List<Treatment>>> entry : d_cycleClasses.entrySet()) {
			InconsistencyParameter param = getInconsistencyParameter(entry.getKey());
			if (param != null) {
				parameters.add(param);
			}
		}
		Collections.sort(parameters, new ParameterComparator());
		return parameters;
	}

	@Override
	protected Map<NetworkParameter, Integer> parameterizeFunctional(Treatment ta, Treatment tb) {
		Map<NetworkParameter, Integer> pmtz = super.parameterizeFunctional(ta, tb);
		
		List<Treatment> cycle = GraphUtil.findPath(d_tree, ta, tb);
		cycle.add(ta);
		Partition partition = findPartition(standardizeCycle(cycle));
		if (partition != null) {
			InconsistencyParameter param = getInconsistencyParameter(partition);
			if (param != null) {
				pmtz.put(param, -relativeDirection(rebaseCycle(cycle), param.getCycle()));
			}
		}
		
		return pmtz;
	}

	private int relativeDirection(List<Treatment> cycle, List<Treatment> ref) {
		Treatment plus = ref.get(1); // if encountered first, cycles are in the same direction
		Treatment minus = ref.get(ref.size() - 2); // in encountered first, cycles are in opposite directions
		
		for (int i = 1; i < cycle.size() - 1; ++i) {
			if (cycle.get(i).equals(plus)) {
				return 1;
			} else if (cycle.get(i).equals(minus)) {
				return -1;
			}
		}

		throw new IllegalStateException("Could not get relative directionality for " + cycle + " and " + ref);
	}

	/**
	 * Get the inconsistency parameter associated with this partition.
	 * @param p The partition.
	 * @return An inconsistency parameter, or null if this partition is not potentially inconsistent.
	 */
	private InconsistencyParameter getInconsistencyParameter(Partition p) {
		if (!isInconsistencyCycle(p)) {
			return null;
		}
		
		return new InconsistencyParameter(standardizeCycle(p.asCycle()));
	}
	
	/**
	 * Find the partition this cycle belongs to.
	 * @param cycle A cycle, standardized.
	 * @return The partition, or null if not found.
	 */
	private Partition findPartition(List<Treatment> cycle) {
		for (Entry<Partition, Set<List<Treatment>>> entry : d_cycleClasses.entrySet()) {
			if (entry.getValue().contains(cycle)) {
				return entry.getKey();
			}
		}
		return null;
	}

}
