import gmsh

# Initialize gmsh
gmsh.initialize()

# Create a new model
gmsh.model.add("rectangle_with_notch_and_hole")

# Define the rectangle dimensions
length = 20
height = 8

# Define the notch dimensions and position
notch_length = 1
notch_position = 4

# Element type
isQuadMesh = True

# Define the circular hole dimensions and position
hole_radius = 0.25
hole_center_x = 6
hole_center_y = 2.75
hole_center_y2 = 2.75 + 2
hole_center_y3 = 2.75 + 4

# Create the main rectangle
rectangle = gmsh.model.occ.addRectangle(0, 0, 0, length, height)
# Fragment the rectangle into three parts
fragment_points = [
    (3, 0, 0),
    (7.5, 0, 0)
]

# Add points to the model
fragment_point_tags = [gmsh.model.occ.addPoint(x, y, z) for x, y, z in fragment_points]

# Create vertical lines for fragmentation
fragment_lines = [
    gmsh.model.occ.addLine(fragment_point_tags[0], gmsh.model.occ.addPoint(3, height, 0)),
    gmsh.model.occ.addLine(fragment_point_tags[1], gmsh.model.occ.addPoint(7.5, height, 0))
]

# Fragment the rectangle with the vertical lines
recfrag = gmsh.model.occ.fragment([(2, rectangle)], [(1, fragment_lines[0]), (1, fragment_lines[1])])
first_rec_tag = recfrag[1][0][0][1]
second_rec_tag = recfrag[1][0][1][1]
third_rec_tag = recfrag[1][0][2][1]
# print (second_rec_tag)
# exit()
gmsh.model.occ.synchronize()

# Create a rectangle starting at (5,0) and ending at (8,2)
rectangle2 = gmsh.model.occ.addRectangle(5.5, 0, 0, 1.9, 1.5)
# Merge this rectangle to second_rec_tag
recfrag2 = gmsh.model.occ.fragment([(2, second_rec_tag)], [(2, rectangle2)])
first_rec2 = recfrag2[1][0][0][1]
second_rec2 = recfrag2[1][0][1][1]

gmsh.model.occ.synchronize()

# Create a rectangle starting at (5,0) and ending at (8,2)
rectangle3 = gmsh.model.occ.addRectangle(5.0, 7.5, 0, 2.4, 0.5)
# Merge this rectangle to second_rec_tag
recfrag3 = gmsh.model.occ.fragment([(2, first_rec2)], [(2, rectangle3)])
first_rec3 = recfrag3[1][0][0][1]
second_rec3 = recfrag3[1][0][1][1]

gmsh.model.occ.synchronize()


# Create the triangular notch
notchsize = 0.025
notch_points = [
    (notch_position - notchsize, 0, 0),  # Left point of the base
    (notch_position + notchsize, 0, 0),  # Right point of the base
    (notch_position, notch_length, 0)        # Top point of the triangle
]

# Add points to the model
notch_point_tags = [gmsh.model.occ.addPoint(x, y, z) for x, y, z in notch_points]

# Create lines for the triangular notch
notch_lines = [
    gmsh.model.occ.addLine(notch_point_tags[0], notch_point_tags[1]),  # Base
    gmsh.model.occ.addLine(notch_point_tags[1], notch_point_tags[2]),  # Right side
    gmsh.model.occ.addLine(notch_point_tags[2], notch_point_tags[0])   # Left side
]

# Create a curve loop and plane surface for the triangular notch
notch_curve_loop = gmsh.model.occ.addCurveLoop(notch_lines)
notch = gmsh.model.occ.addPlaneSurface([notch_curve_loop])

# Create the circular hole
hole = gmsh.model.occ.addDisk(hole_center_x, hole_center_y, 0, hole_radius, hole_radius)
# Create the second circular hole
hole2 = gmsh.model.occ.addDisk(hole_center_x, hole_center_y2, 0, hole_radius, hole_radius)

# Create the third circular hole
hole3 = gmsh.model.occ.addDisk(hole_center_x, hole_center_y3, 0, hole_radius, hole_radius)

# Cut the notch and hole from the rectangle
gmsh.model.occ.cut([(2, first_rec3)], [(2, notch), (2, hole), (2, hole2), (2, hole3)], removeObject=True, removeTool=True)

# Synchronize the CAD kernel with the Gmsh model
gmsh.model.occ.synchronize()

# Define physical groups for the points
bottom_left_point = gmsh.model.occ.addPoint(1, 0, 0)
bottom_right_point = gmsh.model.occ.addPoint(length - 1, 0, 0)
top_middle_point = gmsh.model.occ.addPoint(length / 2, height, 0)

# Fragment the geometry to ensure points are part of the lines and surfaces
gmsh.model.occ.fragment([(0, bottom_left_point)], [(2, first_rec_tag)])
gmsh.model.occ.fragment([(0, bottom_right_point), (0, top_middle_point)], [(2, third_rec_tag)])

# Synchronize the CAD kernel with the Gmsh model
gmsh.model.occ.synchronize()


# Add physical groups for the points
bottom_left_group = gmsh.model.addPhysicalGroup(0, [bottom_left_point])
gmsh.model.setPhysicalName(0, bottom_left_group, "ptfixedxy")

bottom_right_group = gmsh.model.addPhysicalGroup(0, [bottom_right_point])
gmsh.model.setPhysicalName(0, bottom_right_group, "ptfixedy")

top_middle_group = gmsh.model.addPhysicalGroup(0, [top_middle_point])
gmsh.model.setPhysicalName(0, top_middle_group, "ptdispy")

# Synchronize the CAD kernel with the Gmsh model
gmsh.model.occ.synchronize()

# Add a physical group for the 2D domain
domain_group = gmsh.model.addPhysicalGroup(2, [first_rec_tag,third_rec_tag,second_rec2,second_rec3])
gmsh.model.setPhysicalName(2, domain_group, "dom")

domain_group_frac = gmsh.model.addPhysicalGroup(2, [first_rec3])
gmsh.model.setPhysicalName(2, domain_group_frac, "domfrac")

# Synchronize the CAD kernel with the Gmsh model
gmsh.model.occ.synchronize()

# Show the model in the Gmsh graphical user interface
gmsh.fltk.run()

# Define a field for mesh refinement
field_id = gmsh.model.mesh.field.add("Box")
gmsh.model.mesh.field.setNumber(field_id, "VIn", 0.03)  # Target mesh size inside the box
gmsh.model.mesh.field.setNumber(field_id, "VOut", 1.25)   # Target mesh size outside the box
gmsh.model.mesh.field.setNumber(field_id, "XMin", 3)
gmsh.model.mesh.field.setNumber(field_id, "XMax", 7.5)
gmsh.model.mesh.field.setNumber(field_id, "YMin", 0)
gmsh.model.mesh.field.setNumber(field_id, "YMax", height)
gmsh.model.mesh.field.setNumber(field_id, "ZMin", 0)
gmsh.model.mesh.field.setNumber(field_id, "ZMax", 0)

# Set the mesh size field as the background mesh field
gmsh.model.mesh.field.setAsBackgroundMesh(field_id)

# Set the mesh algorithm to use quadrilateral elements
if isQuadMesh == True:
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.RecombineAll", 1)

# Generate the mesh
gmsh.model.mesh.generate(2)

# Show the model in the Gmsh graphical user interface
gmsh.fltk.run()

# Save the mesh to a file
if isQuadMesh == True:
    gmsh.write("bittencourt_quad.msh")
else:
    gmsh.write("bittencourt.msh")

# Finalize gmsh
gmsh.finalize()