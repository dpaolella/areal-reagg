package main

import (
	"bytes"
	"fmt"
	"strconv"
	"strings"

	"github.com/ctessum/geom"
	"github.com/ctessum/geom/encoding/shp"
	"github.com/ctessum/geom/index/rtree"
)

// This is the spatial projection of the InMAP grid.
const gridProj = `PROJCS["Lambert_Conformal_Conic",GEOGCS["GCS_unnamed ellipse",DATUM["D_unknown",SPHEROID["Unknown",6370997,0]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic"],PARAMETER["standard_parallel_1",33],PARAMETER["standard_parallel_2",45],PARAMETER["latitude_of_origin",40],PARAMETER["central_meridian",-97],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]`

type fromGeom struct {
	geom.Polygon
	fields map[string]float64
}

type toGeom struct {
	geom.Polygonal
	fields map[string]float64
}

func loadFromShp(path string) (records []fromGeom) {

	// Load geometries.
	d, err := shp.NewDecoder(path)
	if err != nil {
		panic(err)
	}

	var fieldNames []string
	for _, f := range d.Fields() {
		fieldNames = append(fieldNames, shpFieldName2String(f.Name))
	}

	// Decode a record from the input file.
	for {

		var rec fromGeom
		g, f, more := d.DecodeRowFields(fieldNames...)
		if !more {
			break
		}

		rec.Polygon = g.(geom.Polygon)
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

	return records

}

func loadToShp(path string) (bounds []*geom.Bounds, polygonTree *rtree.Rtree) {

	// Load geometries.
	d, err := shp.NewDecoder(path)
	if err != nil {
		panic(err)
	}

	polygonTree = rtree.NewTree(25, 50) // Create a spatial index for fast searching.
	for {
		var rec toGeom
		if more := d.DecodeRow(&rec); !more {
			break
		}
		polygonTree.Insert(rec) // insert the urban area into the spatial index.
		bounds = append(bounds, rec.Bounds())
	}

	// Check to see if any errors occured during decoding.
	if err := d.Error(); err != nil {
		panic(err)
	}

	return bounds, polygonTree

}

func main() {

	weightVars := make(map[string]string)
	weightCheck := make(map[string]bool)

	weightCheck["TotalPop"] = false
	weightCheck["TotalPM25"] = true
	weightVars["TotalPM25"] = "TotalPop"

	fmt.Println("Loading origin shapefile")
	records := loadFromShp("2005nei_output_1.shp")
	fmt.Println("Loading destination shapefile")
	bounds, polygonTree := loadToShp("states.shp")

	for _, r := range records {
		rArea := r.Area()
		for _, p := range polygonTree.SearchIntersect(r.Bounds()) {
			isect := r.Intersection(p.(toGeom).Polygonal)
			isectArea := isect.Area()
			p.(toGeom).fields = make(map[string]float64, len(r.fields))
			for k, _ := range r.fields {
				if _, ok := weightCheck[k]; ok {
					weightVar := weightVars[k]
					fmt.Println("Hello World")
					p.(toGeom).fields[k] += r.fields[k] * r.fields[weightVar]
				} else {
					fmt.Println("Hello World")
					p.(toGeom).fields[k] += r.fields[k] * (isectArea / rArea)
				}
			}
		}
	}
	for _, b := range bounds {
		for _, p := range polygonTree.SearchIntersect(b) {
			for k, _ := range p.(toGeom).fields {
				if _, ok := weightCheck[k]; ok {
					weightVar := weightVars[k]
					p.(toGeom).fields[k] /= p.(toGeom).fields[weightVar]
				}
			}
		}
	}
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
