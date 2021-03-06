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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import org.drugis.mtc.graph.GraphUtil;
import org.drugis.mtc.model.Measurement;
import org.drugis.mtc.model.Network;
import org.drugis.mtc.model.Study;
import org.drugis.mtc.model.Treatment;
import org.junit.Test;

import edu.uci.ics.jung.algorithms.transformation.FoldingTransformerFixed.FoldedEdge;
import edu.uci.ics.jung.graph.Hypergraph;
import edu.uci.ics.jung.graph.UndirectedGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class NetworkModelTest {
	@Test
	public void testStudyGraphEmpty() {
		Network network = new Network();
		Hypergraph<Treatment, Study> graph = NetworkModel.createStudyGraph(network);
		assertEquals(0, graph.getVertexCount());
		assertEquals(0, graph.getEdgeCount());
	}
	
	@Test
	public void testStudyGraph() {
		Network network = new Network();
		Treatment ta = new Treatment("A", null);
		Treatment tb = new Treatment("B", null);
		Treatment tc = new Treatment("C", null);
		Treatment td = new Treatment("D", null);
		List<Treatment> treatments = Arrays.asList(ta, tb, tc, td);
		network.getTreatments().addAll(treatments);
		Study s1 = new Study("Study 1");
		s1.getMeasurements().add(new Measurement(ta));
		s1.getMeasurements().add(new Measurement(tb));
		Study s2 = new Study("Study 2");
		s2.getMeasurements().add(new Measurement(ta));
		s2.getMeasurements().add(new Measurement(tb));
		Study s3 = new Study("Study 3");
		s3.getMeasurements().add(new Measurement(ta));
		s3.getMeasurements().add(new Measurement(tb));
		s3.getMeasurements().add(new Measurement(tc));
		Study s4 = new Study("Study 4");
		s4.getMeasurements().add(new Measurement(tb));
		s4.getMeasurements().add(new Measurement(tc));
		List<Study> studies = Arrays.asList(s1, s2, s3, s4);
		network.getStudies().addAll(studies);
		
		Hypergraph<Treatment, Study> graph = NetworkModel.createStudyGraph(network);
		assertEquals(new HashSet<Treatment>(treatments), graph.getVertices());
		assertEquals(new HashSet<Study>(studies), graph.getEdges());
		assertEquals(new HashSet<Treatment>(Arrays.asList(ta, tb)), graph.getIncidentVertices(s1));
		assertEquals(new HashSet<Treatment>(Arrays.asList(ta, tb)), graph.getIncidentVertices(s2));
		assertEquals(new HashSet<Treatment>(Arrays.asList(ta, tb, tc)), graph.getIncidentVertices(s3));
		assertEquals(new HashSet<Treatment>(Arrays.asList(tb, tc)), graph.getIncidentVertices(s4));
		
		assertFalse(GraphUtil.isWeaklyConnected(graph));
	}
	
	@Test
	public void testComparisonGraph() {
		Network network = new Network();
		Treatment ta = new Treatment("A", null);
		Treatment tb = new Treatment("B", null);
		Treatment tc = new Treatment("C", null);
		Treatment td = new Treatment("D", null);
		List<Treatment> treatments = Arrays.asList(ta, tb, tc, td);
		network.getTreatments().addAll(treatments);
		Study s1 = new Study("Study 1");
		s1.getMeasurements().add(new Measurement(ta));
		s1.getMeasurements().add(new Measurement(tb));
		Study s2 = new Study("Study 2");
		s2.getMeasurements().add(new Measurement(ta));
		s2.getMeasurements().add(new Measurement(tb));
		Study s3 = new Study("Study 3");
		s3.getMeasurements().add(new Measurement(ta));
		s3.getMeasurements().add(new Measurement(tb));
		s3.getMeasurements().add(new Measurement(tc));
		Study s4 = new Study("Study 4");
		s4.getMeasurements().add(new Measurement(tb));
		s4.getMeasurements().add(new Measurement(tc));
		List<Study> studies = Arrays.asList(s1, s2, s3, s4);
		network.getStudies().addAll(studies);
		
		UndirectedGraph<Treatment, FoldedEdge<Treatment, Study>> graph = NetworkModel.createComparisonGraph(network);
		assertEquals(new HashSet<Treatment>(treatments), new HashSet<Treatment>(graph.getVertices()));
		assertEquals(3, graph.getEdgeCount());
		assertEquals(new HashSet<Study>(Arrays.asList(s1, s2, s3)), new HashSet<Study>(graph.findEdge(ta, tb).getFolded()));
		assertEquals(new HashSet<Study>(Arrays.asList(s3)), new HashSet<Study>(graph.findEdge(ta, tc).getFolded()));
		assertEquals(new HashSet<Study>(Arrays.asList(s3, s4)), new HashSet<Study>(graph.findEdge(tb, tc).getFolded()));
		
		assertEquals(new Pair<Treatment>(ta, tb), graph.getEndpoints(graph.findEdge(ta, tb)));
		
		assertFalse(GraphUtil.isWeaklyConnected(graph));
	}
}
