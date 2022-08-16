# ifndef __tetramesh__
# define __tetramesh__

# include "../interpolation/point.h"
# include <vector>
# include <map>

// ------------------------------------------------------------------------- //
struct tetra
{
	
	// element id
	int id;					
	// nodes index
	int nodes[4];		
	// center
	point center;		
	// volume
	double volume;	
	// th_boundary element
	bool belm;	
	// ext face, global id
	std::vector<int> ext_faces;
	// associated element set
	std::string elm_set;	

};
// ------------------------------------------------------------------------- //
struct tri
{
	
public:

	int id;										// face id
	std::string face_set;			// associated face set
	int nodes[3];							// nodes index
	point normal;							// normal vector
	point center;							// center
	double area;							// area
	int elm;									// index of the associated element
	// stores the id of this face according to the node numbering
	// of the th_boundary element
	// ex nodes = {0,1,2,3} and the face is defined as 
	// faces = {3,2,0}
	// then face_id will be 1
	int elm_fid;

};
// ------------------------------------------------------------------------- //
typedef std::map<std::string, std::vector<int>> geom_set;
// ------------------------------------------------------------------------- //
struct tetra_mesh
{
	
	std::vector<point> nodes;
	std::vector<tri> faces;
	std::vector<tetra> elements;
	
	std::map<int,int> o2n; // original to new node numbering
	std::map<int,int> n2o; // new to original node numbering
	
	geom_set	node_sets;
	geom_set	face_sets;
	geom_set	elm_sets;

	void read_abaqus(const char* filename, double scale_factor);
	void check_free_nodes();
	void define_elements();
	void define_faces();
	void ns2fs();
	void build(const char* filename, double scale);
	void write_fset_vtk(const char* fset) const;


};
// ------------------------------------------------------------------------- //
struct tetra_sub_mesh
{
	const tetra_mesh& msh;

	std::vector<std::string> sets;		// list of volume sets
	std::vector<int> gnid;						// list of global nodes id
	std::vector<int> geid;						// list of global elm id
	std::map<int,int>	g2l_nn;					// global to local node numbering

	int ne() const { return geid.size(); }
	int nn() const { return gnid.size(); }
	const point& get_node(int lid) const { return msh.nodes[gnid[lid]]; }
	const tetra& get_elm(int lid) const { return msh.elements[geid[lid]]; }

	// just initialize the refernece to the original mesh
	tetra_sub_mesh(const tetra_mesh& mesh) : msh(mesh) {}
	// build from a list of volumes set
	void build(std::vector<std::string> set_list);
	// write the subset in vtk format
	void write_vtk(const char* fname);
};

# endif
