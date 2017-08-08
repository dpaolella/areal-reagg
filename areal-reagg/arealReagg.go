package arealReagg

import (
	"bytes"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/ctessum/geom"
	"github.com/ctessum/geom/encoding/shp"
	"github.com/ctessum/geom/index/rtree"
	"github.com/ctessum/geom/proj"
	goshp "github.com/jonas-p/go-shp"
)

// This is the spatial projection of the InMAP grid.
const gridProj = `PROJCS["Lambert_Conformal_Conic",GEOGCS["GCS_unnamed ellipse",DATUM["D_unknown",SPHEROID["Unknown",6370997,0]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",33],PARAMETER["standard_parallel_2",45],PARAMETER["latitude_of_origin",40],PARAMETER["central_meridian",-97],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]`

type inRec struct {
	geom.Polygonal
	fields map[string]float64
}

type outRec struct {
	geom.Polygon
	fields map[string]float64
}

func loadFromShp(path string) (records []inRec, fieldNames []string) {

	// Load geometries.
	d, err := shp.NewDecoder(path)
	if err != nil {
		panic(err)
	}

	for _, f := range d.Fields() {
		fieldNames = append(fieldNames, shpFieldName2String(f.Name))
	}

	// Decode a record from the input file.
	for {

		var rec inRec
		g, f, more := d.DecodeRowFields(fieldNames...)
		if !more {
			break
		}

		rec.Polygonal = g.(geom.Polygonal)
		rec.fields = make(map[string]float64, len(f))
		for k, v := range f {
			rec.fields[k], err = s2f(v)
			if err != nil {
				panic(err)
			}
		}
		records = append(records, rec)
	}

	// Check to see if any errors occured during decoding.
	if err = d.Error(); err != nil {
		panic(err)
	}

	return records, fieldNames

}

func loadToShp(path string) (bounds []*geom.Bounds, polygonTree *rtree.Rtree, records map[geom.Bounds]*outRec) {

	// Load geometries.
	d, err := shp.NewDecoder(path)
	if err != nil {
		panic(err)
	}

	sr, err := d.SR() // Get the spatial projection of the shapefile.
	if err != nil {
		panic(err)
	}
	dst, err := proj.Parse(gridProj) // Parse the InMAP spatial projection.
	if err != nil {
		panic(err)
	}
	t, err := sr.NewTransform(dst) // Create a spatial projection transform.
	if err != nil {
		panic(err)
	}
	fmt.Println(sr)
	fmt.Println(t)

	polygonTree = rtree.NewTree(25, 50) // Create a spatial index for fast searching.
	records = make(map[geom.Bounds]*outRec)
	for {
		var rec outRec
		if more := d.DecodeRow(&rec); !more {
			break
		}
		gg, err := rec.Polygon.Transform(t) // Reproject the polygon to the grid spatial reference.
		if err != nil {
			panic(err)
		}
		rec.Polygon = gg.(geom.Polygon)
		polygonTree.Insert(rec)
		bounds = append(bounds, rec.Polygon.Bounds())
		records[*rec.Bounds()] = &rec
	}

	// Check to see if any errors occured during decoding.
	if err := d.Error(); err != nil {
		panic(err)
	}

	return bounds, polygonTree, records

}

func Run() {

	weightVars := make(map[string]string)
	weightVars["TotalPM25"] = "TotalPop"

	fmt.Println("Loading origin shapefile")
	records, fieldNames := loadFromShp("2005nei_output_1.shp")
	fmt.Println("Loading destination shapefile")
	_, polygonTree, outData := loadToShp("statesNew.shp")

	for _, r := range outData {
		r.fields = make(map[string]float64, len(records[0].fields))
	}

	for _, r := range records {
		// if i > 10 {
		// 	break
		// }
		rArea := r.Area()
		fmt.Println(len(polygonTree.SearchIntersect(r.Bounds())))
		for _, p := range polygonTree.SearchIntersect(r.Bounds()) {
			isect := r.Intersection(p.(outRec).Polygon)
			isectArea := isect.Area()
			// outData[*p.Bounds()].fields = make(map[string]float64, len(r.fields))
			for k, _ := range r.fields {
				// fmt.Printf("string: %v\n", k)
				// fmt.Printf("val: %v\n", v)
				if w, ok := weightVars[k]; ok {
					// fmt.Println(w)
					// fmt.Println(r.fields)
					// fmt.Println(outData[*p.Bounds()].fields[k])
					outData[*p.Bounds()].fields[k] += r.fields[k] * r.fields[w]
				} else {
					outData[*p.Bounds()].fields[k] += r.fields[k] * (isectArea / rArea)
				}
			}
		}
	}
	// fmt.Println("is it still populated?")
	// for i, r := range records {
	// 	fmt.Println(i)
	// 	if i > 30 {
	// 		break
	// 	}
	// 	for _, p := range polygonTree.SearchIntersect(r.Bounds()) {
	// 		fmt.Println(outData[*p.Bounds()].fields)
	// 	}
	// }

	for g, v := range outData {
		for k, _ := range v.fields {
			if w, ok := weightVars[k]; ok {
				outData[g].fields[k] /= outData[g].fields[w]
			}
		}
	}

	// fmt.Println("print bounds:")
	// for _, b := range bounds {
	// 	fmt.Println(len(polygonTree.SearchIntersect(b)))
	// 	for _, p := range polygonTree.SearchIntersect(b) {
	// 		for k, _ := range outData[*p.Bounds()].fields {
	// 			if w, ok := weightVars[k]; ok {
	// 				outData[*p.Bounds()].fields[k] /= p.(outRec).fields[w]
	// 			}
	// 		}
	// 	}
	// }
	var outRecs []*outRec
	for _, g := range outData {
		outRecs = append(outRecs, g)
	}

	fields := make([]goshp.Field, len(fieldNames))
	for i, v := range fieldNames {
		fields[i] = goshp.FloatField(v, 14, 8)
	}
	shape, err := shp.NewEncoderFromFields("test.shp", goshp.POLYGON, fields...)
	if err != nil {
		fmt.Printf("error creating output shapefile: %v", err)
	}
	for i, r := range outRecs {
		outFields := make([]interface{}, len(fieldNames))
		for j, v := range fieldNames {
			outFields[j] = outRecs[i].fields[v]
		}
		err = shape.EncodeFields(r.Polygon, outFields...)
		if err != nil {
			fmt.Printf("error writing output shapefile: %v", err)
		}
	}
	shape.Close()

	// Create .prj file
	f, err := os.Create("test.prj")
	if err != nil {
		fmt.Printf("error creating output prj file: %v", err)
	}
	fmt.Fprint(f, gridProj)
	f.Close()

}

func shpFieldName2String(name [11]byte) string {
	b := bytes.Trim(name[:], "\x00")
	n := bytes.Index(b, []byte{0})
	if n == -1 {
		n = len(b)
	}
	return strings.TrimSpace(string(b[0:n]))
}

func s2f(s string) (float64, error) {
	if removeNull(s) == "************************" { // Null value
		return 0., nil
	}
	f, err := strconv.ParseFloat(removeNull(s), 64)
	return f, err
}

func removeNull(s string) string {
	s = s[0 : len(s)-strings.Count(s, "\x00")]
	return s
}
