package org.drugis

class Graph[T <% Ordered[T]](edges: Set[(T, T)]) {
	val edgeSet: Set[(T, T)] = edges
	val vertexSet: Set[T] = vertices(edges)

	def union(other: Graph[T]) =
		if (other canEqual this) new Graph[T](edgeSet ++ other.edgeSet)
		else throw new IllegalArgumentException
	def intersection(other: Graph[T]) = 
		if (other canEqual this) new Graph[T](edgeSet ** other.edgeSet)
		else throw new IllegalArgumentException

	def edgeVector: List[(T, T)] =
		edgeSet.toList.sort((a, b) =>
			if (a._1 == b._1) a._2 < b._2
			else a._1 < b._1)

	def incidenceVector(edgeVector: List[(T, T)]): List[Boolean] = {
		edgeVector.map(e => containsEdge(e))
	}

	def containsEdge(e: (T, T)): Boolean = {
		edgeSet.contains(e)
	}

	override def equals(other: Any) = other match {
		case that: Graph[_] => 
			(that canEqual this) &&
			this.edgeSet == that.edgeSet
		case _ => false
	}

	def canEqual(other: Any) = other.isInstanceOf[Graph[_]]

	override def hashCode = edgeSet.hashCode

	private def vertices(edges: Set[(T, T)]): Set[T] = {
		edges.flatMap(e => List(e._1, e._2))
	}
}

class UndirectedGraph[T <% Ordered[T]](edges: Set[(T, T)])
extends Graph[T](edges.map(e => UndirectedGraph.order(e))) {

	override def union(other: Graph[T]): UndirectedGraph[T] =
		if ((this canEqual other) && (other canEqual this))
			new UndirectedGraph[T](edgeSet ++ other.edgeSet)
		else throw new IllegalArgumentException

	override def intersection(other: Graph[T]): UndirectedGraph[T] =
		if ((this canEqual other) && (other canEqual this))
			new UndirectedGraph[T](edgeSet ** other.edgeSet)
		else throw new IllegalArgumentException

	override def containsEdge(e: (T, T)): Boolean =
		super.containsEdge(UndirectedGraph.order(e))

	override def equals(other: Any) = other match {
		case that: UndirectedGraph[_] =>
			(that canEqual this) &&
			this.edgeSet == that.edgeSet
		case _ => false
	}

	override def canEqual(other: Any) = other.isInstanceOf[UndirectedGraph[_]]
}

object UndirectedGraph {
	def order[T <% Ordered[T]](e: (T, T)): (T, T) = 
		if (e._1 <= e._2) e
		else (e._2, e._1)
}