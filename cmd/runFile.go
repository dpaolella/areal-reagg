// Copyright © 2017 NAME HERE <EMAIL ADDRESS>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package cmd

import (
	"github.com/dpaolella/areal-reagg/areal-reagg"
	"github.com/spf13/cobra"
)

// runFileCmd represents the runFile command
var runFileCmd = &cobra.Command{
	Use:   "runFile",
	Short: "Run on a single shapefile",
	Long: `runFile exports data from a shapefile to a
	a csv of the same name`,
	Run: func(cmd *cobra.Command, args []string) {
		arealReagg.Run()
	},
}

func init() {
	RootCmd.AddCommand(runFileCmd)
}
