val r1 = List[Boolean](false, true)
val r3 = List[Boolean](true, true)
val r2 = List[Boolean](true, false)
val rows = List[List[Boolean]](r1, r2, r3)

val matrix = new Matrix(rows)

println(matrix)
println(gaussElimGF2(matrix))
println(gaussElimGF2(matrix).elements.findIndexOf(l => !l.contains(true)))