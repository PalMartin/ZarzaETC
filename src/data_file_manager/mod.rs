//
// data_file_manager/main.rs: Handle data files and defult paths
//
// Copyright (c) 2025 Pablo Álvarez Martín <pablo.alvmar12@gmail.com>
// Copyright (c) 2025 Gonzalo J. Carracedo <BatchDrake@gmail.com>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
//

// LIBRARIES IM FINDING NECESSARY TO DO THIS, CHECK IF THEY ARE THE CORRECT OF MORE OPTIMAL ONES
// nix::unitsd -> Safe wrappers around functions found in libc “unistd.h” header
//use nix::unistd::{AccessFlags, access}; -> apparently, these are no longer in 'unistd', now they are in fcntl
//use nix::fcntl::{AccessFlags, access};
use nix::errno::Errno; // errno() from nix::errno::errno is deprecated, best use errno::errno;
use std::ffi::CString; 
//use nix::libc::{access, strdup, W_OK, R_OK, F_OK, X_OK, stat, S_IFDIR};
use nix::libc::*;
use std::path::Path;
use errno::errno;
use once_cell::sync::Lazy; // for the singleton
use std::env;
//use core_foundation::bundle::{CFBundle, CFBundleRef, CFBundleGetMainBundle, CFBundleCopyBundleURL, CFBundleCopyResourceURL};
//use core_foundation::url::{CFURL, CFURLRef, CFURLCopyFileSystemPath, kCFURLPOSIXPathStyle};
//use core_foundation::string::{CFString, CFStringRef, CFStringEncoding, CFStringGetSystemEncoding, CFStringGetCStringPtr};
use core_foundation::bundle::*;
use core_foundation::url::*;
use core_foundation::string::*;
use std::ptr;
use std::ffi::CStr;

#[derive(Clone)]
pub struct DataFileManager {
    paths: Vec<String>,
}

#[cfg(target_os = "macos")]
fn get_resource_bundle_path(rel_path: CFStringRef) -> *const i8 { 

    let mut main_bundle: CFBundleRef = ptr::null_mut();
    let mut dir_url: CFURLRef = ptr::null_mut();
    let mut dir_path: CFStringRef = ptr::null_mut();
    let mut encmethod: CFStringEncoding = 0;
    let mut path: *const i8 = ptr::null();

    main_bundle = unsafe { CFBundleGetMainBundle() };
    if !main_bundle.is_null() {
        dir_url = unsafe { CFBundleCopyResourceURL(
            main_bundle,
            rel_path,
            ptr::null(), /* resourceType */
            ptr::null() /* subDirName */) };

        if dir_path.is_null() {
            return path;
        }

        dir_path = unsafe { CFURLCopyFileSystemPath(dir_url, kCFURLPOSIXPathStyle) };
        if dir_path.is_null() {
            return path;
        }

        encmethod = unsafe { CFStringGetSystemEncoding() };
        path = unsafe { CFStringGetCStringPtr(dir_path, encmethod) };
    }
    return path;
}

impl DataFileManager {

    pub fn new() -> Self {
        let mut data_file_manager = DataFileManager {
            paths: Vec::new(),
        };

        let working_dir = env::current_dir().expect("Failed to retrieve the current working directory");
        data_file_manager.add_search_path(working_dir.display().to_string());

        let extra_path = env::var("ZARZA_ETC_DATA_DIR");
        if extra_path.is_ok() {
            let path = extra_path.unwrap();
            if !path.is_empty() {
                data_file_manager.add_search_path(path);
            }
        }

        #[cfg(target_os = "macos")]
        {
            let name = CString::new("ZarzaETC").expect("CString::new failed");
            let cf_name = unsafe {
                CFStringCreateWithCString(
                    ptr::null(), 
                    name.as_ptr(), 
                    CFStringGetSystemEncoding()
                ) }; // Convert *const i8 to *const __CFString
            let rsrc_path = get_resource_bundle_path(cf_name);
            if !rsrc_path.is_null() {
                let c_str = unsafe { CStr::from_ptr(rsrc_path) };
                let rsrc_str = c_str.to_str().expect("Failed to convert C string to Rust string");
                data_file_manager.clone().add_search_path(format!("{}/data", rsrc_str));
            }
        }
        return data_file_manager;
    }

    // singleton instance of DataFileManager
    pub fn instance() -> &'static DataFileManager {
        static INSTANCE: Lazy<DataFileManager> = Lazy::new(|| DataFileManager::new());
        &INSTANCE
    }

    pub fn add_search_path(&mut self, path: String) -> bool {

        let mut sbuf: stat = unsafe { std::mem::zeroed() }; // Initialize sbuf with zeroed memory -> apparently this is necessary but i dont understand the whole process

        let c_path = CString::new(path.clone()).expect("CString::new failed");

        if unsafe { stat(c_path.as_ptr(), &mut sbuf as *mut stat) } == -1 {
            let err = errno();
            let err_str = err.to_string();
            eprintln!("warning: cannot stat `{}`: {}", path, err_str);
            return false;
        }

        if sbuf.st_mode & S_IFDIR == 0 {
            eprintln!("warning: file `{}` is not a directory", path);
            return false;
        }

        self.paths.insert(0, path.to_string());
        return true;
    }   

    pub fn find(&self, path: String, flags: i32) -> String {
        if path.is_empty() {
            return path;
        }

        if path.chars().next().unwrap() == '/' {
            let c_path = CString::new(path.clone()).expect("CString::new failed");
            if unsafe { access(c_path.as_ptr(), flags) != -1 } {
                return path;
            }

            // Write access: check if we can create a file there
            if (flags & W_OK) != 0 && Errno::last() == Errno::ENOENT {
                // Copy of path (is what strdup() is doing)
                let copy = CString::new(path.clone()).expect("Failed to create CString"); // check if this is neccesary, i dont think i need to manage the memory manually...s
                if copy.as_ptr().is_null() {
                    panic!("Memory exhausted"); // Equivalent to throwing a runtime error
                }
                let dir = Path::new(&path).parent();
                let dir_str = dir.unwrap().to_str().unwrap();
                let c_dir_str = CString::new(dir_str).expect("Failed to create CString");
                let ret = unsafe { access(c_dir_str.as_ptr(), W_OK | X_OK) };

                if ret != -1 {
                    return path;
                }
            }

        } else {
            for p in &self.paths {
                let full_path = p.clone() + "/" + &path;
                let c_full_path = CString::new(full_path.clone()).expect("CString::new failed");
                if unsafe { access(c_full_path.as_ptr(), flags) } != -1 {
                    return full_path;
                }

                // Write access:: check if we can create a file here
                if (flags & W_OK) != 0 && Errno::last() == Errno::ENOENT {
                    let c_p = CString::new(p.clone()).expect("CString::new failed");
                    if unsafe { access(c_p.as_ptr(), W_OK | X_OK) != -1 } {
                        return full_path;
                    }
                }
            }
        }
        // If nothing passes, return the original path.
        return path;
    }

    pub fn resolve(&self, path: String) -> String {
        return self.find(path, R_OK)
    }

    pub fn suggest(&self, path: String) -> String {
        return self.find(path, W_OK)
    }

    pub fn search_paths(&self) -> Vec<String> {
        return self.paths.clone();
    }

    pub fn data_file(path: &String) -> String {

        let result = DataFileManager::instance().resolve(path.clone());
        if result.is_empty() {
            eprintln!("Required datafile. {} not found.", path);
        }

        return result;
    }
    
}

